

#some useful functions modified from Mudpy
def build_TauPyModel(out_directory,vel_mod_file,background_model='PREM'):
    '''
    This function will take the structure from the .mod file
    and paste it on top of a pre computed mantle structure such as PREM.
    This assumes that the .mod file provided by the user ends with a 0 thickness 
    layer on the MANTLE side of the Moho
    '''
    from numpy import genfromtxt
    from os import environ,path
    from obspy.taup import taup_create, TauPyModel
    import obspy
    #mudpy source folder

    #load user specified .mod infromation
    structure = genfromtxt(vel_mod_file)

    #load background velocity structure
    if background_model=='PREM':

        bg_model_file=obspy.__path__[0]+'/taup/data/prem.nd'
        #Q values
        Qkappa=1300
        Qmu=600

        #Write new _nd file one line at a time
        nd_name=path.basename(vel_mod_file).split('.')[0]
        nd_name=nd_name+'.nd'
#        f=open(out_directory+'/structure/'+nd_name,'w')
        f=open(out_directory+'/'+nd_name,'w')
        #initalize
        ztop=0
        for k in range(len(structure)-1):
            #Variables for writing to file
            zbot=ztop+structure[k,0]
            vp=structure[k,2]
            vs=structure[k,1]
            rho=structure[k,3]
            # Write to the file
            line1=('%8.2f\t%8.5f   %7.5f   %7.5f\t%6.1f     %5.1f\n' % (ztop,vp,vs,rho,Qkappa,Qmu))
            line2=('%8.2f\t%8.5f   %7.5f   %7.5f\t%6.1f     %5.1f\n' % (zbot,vp,vs,rho,Qkappa,Qmu))
            f.write(line1)
            f.write(line2)
            #update
            ztop=zbot
        #now read PREM file libe by libne and find appropriate depth tos tart isnerting
        fprem=open(bg_model_file,'r')
        found_depth=False
        while True:
            line=fprem.readline()
            if line=='': #End of file
                break
            if found_depth==False:
                #Check that it's not a keyword line like 'mantle'
                if len(line.split())>1:
                    #not a keyword, waht depth are we at?
                    current_depth=float(line.split()[0])
                    if current_depth > zbot: #PREM depth alrger than .mod file last line
                        found_depth=True
                        f.write('mantle\n')
            #Ok you have found the depth write it to file
            if found_depth == True:
                f.write(line)
        fprem.close()
        f.close()
        # make TauPy npz
#        taup_in=home+project_name+'/structure/'+nd_name
#        taup_out=home+project_name+'/structure/'
        taup_in=out_directory+'/'+nd_name
        taup_out=out_directory+'/'
        taup_create.build_taup_model(taup_in,output_folder=taup_out)
    else: #To be done later (ha)
        print('ERROR: That background velocity model does not exist')

