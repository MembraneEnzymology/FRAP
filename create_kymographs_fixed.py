import os

from loci.formats import ImageReader
from ij import IJ
from ij.gui import Line
from java.util import Properties
from java.lang import System
from ij.plugin import ChannelSplitter


path = 'G:\\My Drive\\EXPERIMENT DATA\\CONFOCAL LSM 710\\2021\\20210713\\M smegmatis mc2 155 pSMT3 dsRed\\FRAP16X16' 
profile_width = 3


# globally sets the line width of ImageJ
Line.setWidth(profile_width)


def open_image(path):
	"""
	This function opens an image file. If it is an image stack we open the first channel that is not T PMT
	"""
	imp = IJ.openImage(path)
	
	if imp.isHyperStack():
		imps = ChannelSplitter.split(imp)
		for imp in imps:
			label = imp.getStack().getSliceLabel(1)
			if label != 'T PMT':
				return imp
		return imps[0]
	else:
		return imp
		

def read_metadata(path):
	reader = ImageReader()
	reader.setId(path)
	channels = reader.getEffectiveSizeC()
	count = reader.getImageCount() / channels
	metadata = reader.getSeriesMetadata()
	timestamps = [0] * count
	
	for key, value in metadata.items():
		if key.startswith('TimeStamp'):
			position = int(key[11:]) - 1
			timestamps[position] = float(value)
		elif key.startswith('VoxelSizeX'):
			pixel_size = float(value)


	return timestamps, pixel_size


def read_selections(path):
	selections = list()
	
	with open(path) as f:
		lines = f.readlines()
		columns = [column.strip() for column in lines[0].split(',')]
		
		for line in lines[1:]:
			try:
				values = dict(zip(columns, [value.strip() for value in line.split(',')]))
				selections.append((
					values['filename'],
					int(values['photo_bleached']),
					float(values['background']),
					float(values['x1']),
					float(values['y1']),
					float(values['x2']),
					float(values['y2']),
					float(values['aspect_ratio'])
				))
			except:
				pass
			
	return selections


def get_profiles(path, selection):
	filename, photo_bleached, background, x1, y1, x2, y2, aspect_ratio = selection
	image_path = os.path.join(path, filename)

	if not os.path.exists(image_path):
		print('image \'%s\' does not exist, image will be skipped.')
		return

	imp = open_image(image_path)
	n = imp.getImageStackSize();
	reference = list()
	profiles = list()
	line = Line(x1, y1, x2, y2)
	
	for i in range(1, photo_bleached):
		imp.setSlice(i)
		imp.setRoi(line)
		profile = line.getPixels()
		
		if not reference:
			reference = profile
		else:
			reference = [a + b for a, b in zip(reference, profile)]

	# average profile of all slice before photo bleaching
	reference = [float(a) / (photo_bleached - 1) for a in reference]
	
	# and subtract backgound
	reference = [a - background for a in reference]

	# determine the integrated value (sum) of reference
	reference_integrated = sum(reference)
	
	for i in range(photo_bleached, n + 1):
		imp.setSlice(i)
		imp.setRoi(line)

		profile = line.getPixels()
		
		# subtract background
		profile = [a - background for a in profile]

		# determine integrated value
		profile_integrated = sum(profile)

		# divide by integrated profile
		# multiply with integrated reference
		# divide by reference profile
		profile = [a / profile_integrated for a in profile]
		profile = [a * reference_integrated for a in profile]
		profile = [a / b for a, b in zip(profile, reference)]
		
		profiles.append(profile)

	# read timestamps and pixel size
	timestamps, pixel_size = read_metadata(image_path)

	# remove all the time stamps before photo bleaching and make sure they start from 0
	timestamps = timestamps[photo_bleached:]
	timestamps = [t - timestamps[0] for t in timestamps]
	
	return profiles, timestamps, pixel_size


def write_kymograph(path, profiles, timestamps, pixel_size):
	with open(path, 'w') as f:
		distance = ['%f' % (i * pixel_size) for i in range(len(profiles[0]))]
		str = 'time,' + ','.join(distance) + '\n'
		for timestamp, profile in zip(timestamps, profiles):
			str += '%f,' % (timestamp)
			str += ','.join(['%f' % value for value in profile]) + '\n'
		f.write(str)
	
	
def process(path):
	selections_path = os.path.join(path, 'selections.csv')

	if not os.path.exists(selections_path):
		print('selections.csv could not be found in \'%s\'. (re-)run create_selections.py' % path)
		return
	
	selections = read_selections(selections_path)

	for selection in selections:
		print('processing %s' % selection[0])
		profiles, timestamps, pixel_size = get_profiles(path, selection)
		kymograph_path = os.path.join(path, selection[0][:-4] + '-values.csv')
		write_kymograph(kymograph_path, profiles, timestamps, pixel_size)

process(path)



