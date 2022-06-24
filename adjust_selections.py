import os
import math

from ij import IJ
from ij.plugin import ZProjector
from ij.plugin import Zoom
from ij.plugin.frame import RoiManager
from ij.gui import Line
from ij.gui import WaitForUserDialog
from ij.gui import YesNoCancelDialog
from ij.plugin import ChannelSplitter




path = 'G:\\My Drive\\EXPERIMENT DATA\\CONFOCAL LSM 710\\2021\\20210713\\M smegmatis mc2 155 pSMT3 dsRed\\FRAP16X16\\' 
profile_width = 3

# globally sets the line width of ImageJ
Line.setWidth(profile_width)


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
	
	selections_path = os.path.join(path, 'selections.csv')
	
	if not os.path.exists(selections_path):
		print('run create_selections.py to create selections.csv file')
		return

	selections = read_selections(selections_path)
	trimmed_selections = list()
	
	images = list()
	
	for selection in selections:
		image_path = os.path.join(path, selection[0])

		if not os.path.exists(image_path):
			print('skipping %s, file doesn\'t exist' % (selection[0]))
		else:		
			images.append(read_image(image_path, selection[1]))
			trimmed_selections.append(selection)
			
	offsets = create_grid(images)
	set_selections(offsets, trimmed_selections)

	zoom = Zoom()
	for i in range(5):
		zoom.run('in')

	dialog = WaitForUserDialog('Adjust selections and click OK when ready')
	dialog.show()

	dialog = YesNoCancelDialog(None, 'Save selections', 'Do you want to save the selections?')
	
	if dialog.yesPressed():
		adjusted_selections = get_adjusted_selections(trimmed_selections, offsets)
		write_selections(os.path.join(path, 'selections.csv'), adjusted_selections)

	

def read_selections(path):
	selections = list()
	
	with open(path) as f:
		lines = f.readlines()
		
		for line in lines[1:]:
			values = line.split(',')
			selections.append((
				values[0].strip(),	# filename
				int(values[1]),		# photo bleached slice
				float(values[2]),	# background
				float(values[3]),	# x1
				float(values[4]),	# y1
				float(values[5]),	# x2
				float(values[6]),	# y2
				float(values[7])))	# aspect ratio

	return selections


def read_image(path, photo_bleached_slice):
	# open image
	imp = open_image(path)
	
	# create reference
	z_projector = ZProjector(imp)
	z_projector.setMethod(ZProjector.AVG_METHOD)
	z_projector.setStopSlice(photo_bleached_slice - 1)
	z_projector.doProjection()
	imp = z_projector.getProjection()
	
	# normalize
	imp.deleteRoi()
	IJ.run(imp, '32-bit', '')
	statistics = imp.getStatistics()
	ip = imp.getProcessor()

	for y in range(ip.getHeight()):
		for x in range(ip.getWidth()):
			pixel = ip.getf(x, y)
			pixel -= statistics.min
			pixel /= (statistics.max - statistics.min)
			ip.setf(x, y, pixel)
			
	return imp


def create_grid(images):
	# determine number of columns
	columns = int(math.ceil(math.sqrt(len(images))))
	
	image_width = max([imp.getWidth() for imp in images])
	image_height = max([imp.getHeight() for imp in images])
	grid_width = columns * image_width
	grid_height = columns * image_height

	
	# create grid and return (x,y) tuples of each image
	grid_imp = IJ.createImage("selections", grid_width, grid_height, 1, 32)
	grid_ip = grid_imp.getProcessor()
	
	offsets = list()
	
	for i, imp in zip(range(len(images)), images):
		ip = imp.getProcessor()
				
		x0 = (i % columns) * image_width
		y0 = (i // columns) * image_height
		offsets.append((x0, y0))
		
		for y in range(ip.getHeight()):
			for x in range(ip.getWidth()):
				pixel = ip.getf(x, y)
				grid_ip.setf(x0 + x, y0 + y, pixel)

	grid_imp.show()
	
	return offsets


def set_selections(offsets, selections):
	roiManager = RoiManager.getInstance()
	if not roiManager:
		roiManager = RoiManager()

	roiManager.reset()

	for offset, selection in zip(offsets, selections):
		x0, y0 = offset[0], offset[1]
		x1, y1, x2, y2 = selection[3:7]
		line = Line(x0 + x1, y0 + y1, x0 + x2, y0 + y2)
		line.setName(selection[0])
		roiManager.addRoi(line)
		
	roiManager.runCommand('UseNames', 'true')
	roiManager.runCommand('Show All with labels')


def get_adjusted_selections(selections, offsets):
	roiManager = RoiManager.getInstance()

	rois = roiManager.getRoisAsArray()
	adjusted_selections = list()
	
	for selection, offset, roi in zip(selections, offsets, rois):
		x0, y0 = offset
		x1, y1, x2, y2 = roi.x1d - x0, roi.y1d - y0, roi.x2d - x0, roi.y2d - y0
		
		adjusted_selections.append((
			selection[0],
			selection[1],
			selection[2],
			x1,
			y1,
			x2,
			y2,
			selection[7]))
			
	return adjusted_selections
	

def write_selections(selections_path, selections):
	with open(selections_path, 'w') as f:
		f.write('%20s, %10s, %10s, %10s, %10s, %10s, %10s, %10s\n' % ('filename', 'photo_bleached', 'background', 'x1', 'y1', 'x2', 'y2', 'aspect_ratio'))
		for selection in selections:
			f.write('%20s, %10d, %10.3f, %10.3f, %10.3f, %10.3f, %10.3f, %10.3f\n' % selection)
			

process(path)