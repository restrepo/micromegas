import re

def read_micrOMEGAS_output(mO_file):
    dict={}
    file=open(mO_file,'r')
    f=file.readlines()
    for line in f:
        #REGEX ^:begin \d: digit, \s: space, \n: new line, .*: any set, (): group
        br=re.search(r'(^\d\.\d+[Ee][+\-]\d+)\s(.*)\n',line) #branchings
        if br:
            #groups: 2:(.*)-> decay name; 1: (^\d\.\d+[Ee][+\-]\d+)-> value
            dict[br.group(2)]=float(br.group(1))

        tw=re.search(r'(^H\->2\*x).*(\d\.\d+[Ee][+\-]\d+).*\n',line) #totalwidth
        if tw:
            dict[tw.group(1)]=float(tw.group(2))

        DMRD=re.search(r'(Omega)(.*)(\d\.\d+[Ee][+/-]\d+)\n',line)
        if DMRD: dict[DMRD.group(1)]=float(DMRD.group(3))
        DDCS1=re.search(r'(^\s+proton)(\s+SI\s)(\d\.\d+[Ee][+/-]\d+)(.*)',line)
        DDCS2=re.search(r'(^\s+neutron)(\s+SI\s)(\d\.\d+[Ee][+/-]\d+)(.*)',line)
        if DDCS1: dict['DDCS->proton(pb)']=float(DDCS1.group(3))
        if DDCS2: dict['DDCS->neutron(pb)']=float(DDCS2.group(3))

    return dict

if __name__ == '__main__':
    '''micrOMEGAS output in kk'''
    All_data=read_micrOMEGAS_output('kk')
    totalwidth=All_data['H->2*x']
    Relic_Density=All_data['Omega']
    Direct_detection=All_data['DDCS->proton(pb)']
    print totalwidth,All_data['H -> s,S']Relic_Density,Direct_detection
