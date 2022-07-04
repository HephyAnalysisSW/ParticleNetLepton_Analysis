import os
from shutil import copyfile

src='/groups/hephy/cms/suman.chatterjee/Store/Analysis_output/SingleLep/v20/Trees/'
#samples=['SingleMuon','EGamma','TT','singleTop','W_NLO_selb','W_NLO_selc','W_NLO_selj','DY_b_LO','DY_c_LO','DY_l_LO','WW','WZ','ZZ']
samples=['WW']
for sample in samples:
	directory=src+sample
	list_file = os.listdir(directory)
	for file in list_file:
		xx =  file.split('_')
		#yy = '%s'%(xx[1])+'%s'%(xx[2])
		print xx
		yy = '%s'%(xx[1])+'%s'%(xx[2])#+'%s'%(xx[3])
		filename = os.path.join(directory,file)
		new_filename = os.path.join(directory,yy)
		#print filename, new_filename
		copyfile(filename, new_filename)
