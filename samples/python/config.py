'''
Some definitions.
'''
import os

## This definition might not work for everybody. Update with whatever you like.
if 'SAMPLES_DBDIR' in os.environ:
    dbDir = os.environ['SAMPLES_DBDIR'] 
elif os.environ['USER'] in ['dietrich.liko']:
    dbDir = '/users/dietrich.liko/caches/Samples'
elif os.environ['USER'] in ['llechner']:
    dbDir = '/afs/hephy.at/data/llechner01/TTGammaEFT/cache/samples/'
elif os.environ['USER'] in ['lukas.lechner']: #CLIP
    dbDir = '/mnt/hephy/cms/lukas.lechner/TTGammaEFT/cache/samples/'
elif os.environ['USER'] in ['dspitzba', 'dspitzbart']:
    dbDir = '/afs/hephy.at/data/dspitzbart01/nanoAOD/'
elif os.environ['USER'] in ['phussain']:
    dbDir = '/afs/hephy.at/data/cms03/phussain/stopsDilepton/cache/samples/'
elif os.environ['USER'] in ['priya.hussain']:
    dbDir = '/mnt/hephy/cms/priya.hussain/StopsCompressed/cache/samples/'
elif os.environ['USER'] in ['schoef']:
    dbDir = '/afs/hephy.at/data/rschoefbeck01/nanoAOD/'
elif os.environ['USER'] in ['robert.schoefbeck']:
    dbDir = '/users/robert.schoefbeck/caches/Samples'
elif os.environ['USER'] in ['erich.birngruber']:
    dbDir = '/users/robert.schoefbeck/caches/Samples'
elif os.environ['USER'] in ['janik.andrejkovic']:
    dbDir = '/users/janik.andrejkovic/caches/Samples'
elif os.environ['USER'] in ['mdoppler']:
    dbDir = '/afs/hephy.at/data/rschoefbeck01/nanoAOD/'
elif os.environ['USER'] in ['mzarucki']:
    dbDir = '/afs/hephy.at/data/mzarucki02/nanoAOD/caches'
elif os.environ['USER'] in ['suman.chatterjee']:
    dbDir = '/users/suman.chatterjee/caches/Samples'
else:
    dbDir = '/afs/hephy.at/data/%s01/nanoAOD/'%os.environ['USER']

redirector_BE         = 'root://maite.iihe.ac.be/'
redirector_global     = 'root://cms-xrd-global.cern.ch/'
#redirector_clip_local = '/mnt/hephy/cms/'
#redirector_clip_local = '/scratch-cbe/users/hephy/'
#redirector_clip       = 'root://hephyse.oeaw.ac.at:11001/'
redirector_clip       = 'root://eos.grid.vbc.ac.at///eos/vbc/incoming/'
#redirector_clip       = '/eos/vbc/incoming/'

#redirector            = 'root://hephyse.oeaw.ac.at:11001/'
redirector            =  redirector_global 

if not os.path.isdir(dbDir): os.makedirs(dbDir)
