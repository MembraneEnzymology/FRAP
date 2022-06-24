import os

from ij import IJ
from ij.plugin import ZProjector
from ij.plugin import ChannelSplitter



path = 'G:\\My Drive\\EXPERIMENT DATA\\CONFOCAL LSM 710\\2021\\20210713\\M smegmatis mc2 155 pSMT3 dsRed\\FRAP16X16\\' 
bg_selection_width = 5	# for background detection a rectangular selection of this width (and height) is used
selection_line_fraction = 0.8	# change this value to extend or shorten the found selections


def find_photo_bleached_slice(imp):
	"""
	This functions tries to find the first photo bleached frame.
	When the mean intensity of the frame drops below 80% of the mean intensity of the first frame
	it is considered photo bleached. Only the first 10 frames are considered.
	Returns the frame number of the frame that is photo bleached or -1 if no such frame is found.
	"""
	
	stack = imp.getImageStack()
	ip = stack.getProcessor(1)
	statistics = ip.getStatistics()
	threshold = statistics.mean * 0.8
	
	for s in range(2, 11):
		ip = stack.getProcessor(s)
		statistics = ip.getStatistics()
		
		if statistics.mean < threshold:
			return s
	
	return -1


def create_reference(imp, photo_bleached_slice):
	"""
	Averages all frames before photo bleaching.
	"""
	z_projector = ZProjector(imp)
	z_projector.setMethod(ZProjector.AVG_METHOD)
	z_projector.setStopSlice(photo_bleached_slice - 1)
	z_projector.doProjection()
	return z_projector.getProjection()


def measure_background(imp):
	"""
	Measures the background value by taking the corner with the lowest mean pixel intensity.
	"""
	global bg_selection_width
		
	imp.setRoi(0, 0, bg_selection_width, bg_selection_width)
	background = imp.getStatistics().mean
	imp.setRoi(imp.getWidth() - bg_selection_width, 0, bg_selection_width, bg_selection_width)
	background = min(imp.getStatistics().mean, background)
	imp.setRoi(0, imp.getHeight() - bg_selection_width, bg_selection_width, bg_selection_width)
	background = min(imp.getStatistics().mean, background)
	imp.setRoi(imp.getWidth() - bg_selection_width, imp.getHeight() - bg_selection_width, bg_selection_width, bg_selection_width)
	return min(imp.getStatistics().mean, background)


def find_selection(imp):
	"""
	Finds the selection through the cell by fitting an ellipse to the thresholded image.
	"""
	
	imp.deleteRoi()
	imp = imp.duplicate()
	IJ.setAutoThreshold(imp, 'Default dark')
	IJ.run(imp, 'Convert to Mask', '')
	IJ.run(imp, 'Create Selection', '')
	IJ.run(imp, 'Fit Ellipse', '')
	ellipse = imp.getRoi()
	return ellipse.getParams()



def open_image(path):
	"""
	This function opens an image file. If it is an image stack we open the first channel that is not T PMT
	"""
	imp = IJ.openImage(path)
	imp.getTitle()
	
	if imp.isHyperStack():
		imps = ChannelSplitter.split(imp)
		for imp in imps:
			label = imp.getStack().getSliceLabel(1)
			if label != 'T PMT':
				return imp
		return imps[0]
	else:
		return imp

			
def process(path):
	global selection_line_fraction

	filename = os.path.basename(path)

	print('processing %s' % filename)
	
	imp = open_image(path)
	
	photo_bleached_slice = find_photo_bleached_slice(imp)
	
	if photo_bleached_slice == -1:
		print('could not detect photo bleached slice')
		return
	
	reference_imp = create_reference(imp, photo_bleached_slice)
	background = measure_background(reference_imp)
	x1, y1, x2, y2, aspectRatio = find_selection(reference_imp)

	# resize selection
	x0 = (x1 + x2) / 2
	y0 = (y1 + y2) / 2
	dx = ((x2 - x1) / 2) * selection_line_fraction
	dy = ((y2 - y1) / 2) * selection_line_fraction
	x1 = x0 - dx
	y1 = y0 - dy
	x2 = x0 + dx
	y2 = y0 + dy
	
	return filename, photo_bleached_slice, background, x1, y1, x2, y2, aspectRatio


def process_all(path):
	cells = list()

	for filename in os.listdir(path):
		if filename.endswith('.lsm'):
			cells.append(process(os.path.join(path, filename)))

	selections_path = os.path.join(path, 'selections.csv')
	with open(selections_path, 'w') as f:
		f.write('%20s, %10s, %10s, %10s, %10s, %10s, %10s, %10s\n' % ('filename', 'photo_bleached', 'background', 'x1', 'y1', 'x2', 'y2', 'aspect_ratio'))
		for cell in cells:
			f.write('%20s, %10d, %10.3f, %10.3f, %10.3f, %10.3f, %10.3f, %10.3f\n' % cell)

	print(selections_path)

process_all(path)
