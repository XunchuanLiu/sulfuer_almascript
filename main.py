import sys
if '.' not in sys.path:
    sys.path.append('.')
import numpy as np
import os
import glob
import shutil
import astropy.table as t
import casashell

from casatasks import *
#from tclean_cli import tclean_cli as tclean
import casatools.ms
ms = casatools.ms()
tb = casatools.table()

spws = {45: 251.82577, #SO
       }

def run(field='Oph_IRS63', #
        spwid=45,
        doconcatfile=False,
        docubeclean = False,
        scriptconfigfile='scriptconfig.csv',
        mpiscriptdir = './mpicasarun/',
        DELETE=False,
        cubethreshold='10mJy',
        ncore=20,
        mindf = 0.123 #MHz
        ):
    #print("debug")
    sTable = t.Table.read(scriptconfigfile)
    assert field in sTable['source']
    sdex = np.argmax(sTable['source'] == field)
    rawdatadirs = {
        "TM1" : sTable['TM1dir'][sdex],
        "TM2" : sTable['TM2dir'][sdex]
    }
    done = sTable['Done'][sdex]
    
    
    if not os.path.exists('data'):
        os.system('mkdir data')
                
    if os.path.exists('data/%s' %field):
        if DELETE:  #be careful
            if done:
                print("This source has been cleaned!")
                print("If you indeed want to reclean it, set 'done' as 0 for %s in the .csv configure file" %field)
                return 
            else:
                os.system('rm -r data/%s' %field)
                os.system('mkdir data/%s' %field)
    else:
        os.system('mkdir data/%s' %field)
        
    if not os.path.exists(mpiscriptdir):
        os.system('mkdir %s' %mpiscriptdir)
        
    if not (doconcatfile | docubeclean ):
        print('No execution for field %s spwid %i' %(field, spwid))
        return
        
        
    dataroot = '../../archive/2022.1.01411.S/*/*/%s/calibrated'
    if True and doconcatfile:
        splitoutfiles = []
        for configuration in  rawdatadirs:
            rawdatadir = rawdatadirs[configuration]
            memberrootdir = glob.glob(dataroot %rawdatadir)[0] 
            files = glob.glob(os.path.join(memberrootdir,'*_line.ms'))
            print('find %i files for configuration %s' %(len(files),configuration))
            for ifile in range(len(files)):
                assert tb.open(files[ifile])
                splitrawdatacolumn = 'data'
                if 'CORRECTED_DATA' in tb.colnames():
                    splitrawdatacolumn = 'corrected'                
                tb.close()
                theoutputvis=os.path.join("data",field,"%s_spw%i_%s_rawsplit_%i.ms" %(field,spwid,configuration,ifile+1))
                split(vis=files[ifile],
                      spw='%i' %spwid,
                      field=field,
                      outputvis=theoutputvis,
                      datacolumn=splitrawdatacolumn
                      )
                 
                splitoutfiles.append(theoutputvis)      
               

                                
    cleanscripttemp = """tclean(vis='{vis}',
        restart = True,
        calcres = {calcres},
        calcpsf = {calcpsf},
        field = '{field}',
        imagename='{imagename}',
        datacolumn='corrected',
        spw='{spw}',
        specmode='cube', 
        width='{width}',
        start='{start}',
        outframe='LSRK',
        nchan={nchan},
        threshold='{threshold}', 
        imsize=[1200, 1200], 
        cell=['0.036arcsec'], 
        niter={niter}, 
        deconvolver='multiscale',
        scales =[0,5,15,50,150], 
        gridder='mosaic', 
        weighting='briggs',
        robust=0.5,
        pbcor=False, 
        pblimit=0.2, 
        #restoringbeam='common',
        usemask = 'auto-multithresh',
        sidelobethreshold = 2.0,
        noisethreshold = 4.25,
        minbeamfrac = 0.3,
        lownoisethreshold = 1.5, 
        #chanchunks=-1, 
        perchanweightdensity = True,
        interactive=False,
        parallel=True,
        cycleniter={cycleniter},
        nmajor={nmajor},
        negativethreshold=6,
        )
    """
    if True and docubeclean:
        splitoutfiles = glob.glob(os.path.join("data",field,"%s_spw%i_*_rawsplit_*.ms" %(field,spwid)))
        concatfilepath = os.path.join("data",field,"%s_spw%i_rawconcat.ms" %(field,spwid))                       
        if not os.path.exists(concatfilepath):
              print('concat %i .ms files' %len(splitoutfiles))
              concat(vis=splitoutfiles,
                       concatvis=concatfilepath)   

        
        df =None
        fmin = None; fmax = None
        for msfile in splitoutfiles:
            assert ms.open(msfile)
            thef = ms.cvelfreqs(spwids=[0], mode='channel', 
                             width=0, outframe='LSRK')/1E6
            if df is None:
                df = np.abs(thef[1]-thef[0])
            else:
                df = min(  np.abs(thef[1]-thef[0]),df )
            if fmin is None:
                fmin = thef.min()
            else:
                fmin = max(fmin, thef.min()) 
            if fmax is None:
                fmax = thef.max()
            else:
                fmax = min(fmax, thef.max())
            ms.close()
            
            df = max(df,mindf)
        thewidth = '%fMHz' %(df)                 
        thenchan = (fmax-fmin)//df   
        
        theimagename = os.path.join("data",field,'%s_spw%i_cube' %(field,spwid))
        cleanscript = cleanscripttemp.format(
                vis=os.path.join('../',concatfilepath),
                calcres = 'True',
                calcpsf = 'True',
                field = field,
                imagename = os.path.join('../',theimagename),
                spw="",
                width=thewidth,
                start='%fMHz' %fmin,
                nchan = '%i' %thenchan,
                threshold = cubethreshold,
                niter='50000000',
                cycleniter='5000',
                nmajor='5',
            )         
         
        

        scriptfile = os.path.join(mpiscriptdir,'mpiscript_%s_spw%i' %(field,spwid))
        with open(scriptfile,'w') as thefile:
                thefile.write(cleanscript)

    if True and docubeclean:
        #mpiscriptrun=os.path.join(mpiscriptdir,'mpiscriptrun.sh')
        #with open(mpiscriptrun,'w') as thef:
        #    for i in wideband_outids:
        #        scriptfile = 'mpiscript_spw%s' %i
        #        thef.write("mpicasa -n {ncore} casa --nogui --nologger -c {script}\n".format(script=scriptfile, ncore=ncore))
        #open a terminal to execute the mpiscripts, e.g.
        #os.system('gnome-terminal --working-directory="%s" -- bash -c "cat mpiscriptrun.sh"' %(os.path.abspath(mpiscriptdir), ) )
        #However, be aware that gnome-terminal --disable-factory may do not support to open a blocked new terminal
        #so, I use xfce4-terminal instead. It assumes that xfce4-terminal has been installed
        #add --hold if you want to keep the terminal around after the executed command finishes executing
        scriptfile = 'mpiscript_%s_spw%i' %(field,spwid)
        command = "env -u PYTHONPATH -u LD_LIBRARY_PATH mpicasa  -n {ncore} casapy-654p --nogui --nologger  -c {script}  \n".format(script=scriptfile, ncore=ncore) #
        os.system('xfce4-terminal --working-directory="%s" --command "%s"' %(os.path.abspath(mpiscriptdir), command) )

    if True and docubeclean:
            theimagename = os.path.join("data",field,'%s_spw%i_cube.image' %(field,spwid))
            thepb = os.path.join("data",field,'%s_spw%i_cube.pb' %(field,spwid))
            thepbcor = os.path.join("data",field,'%s_spw%i_cube.image.pbcor' %(field,spwid))
            thepbcorfits = os.path.join("data",field,'%s_spw%i_cube.image.pbcor.fits' %(field,spwid))
            impbcor(imagename=theimagename,
                    pbimage=thepb,
                    outfile=thepbcor,
                    overwrite=True,
                   )
            exportfits(imagename=thepbcor,
                       fitsimage=thepbcorfits,
                       overwrite=True, 
                      )    
            
        
    
