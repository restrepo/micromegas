import re
def read_micrOMEGAS_output(mO_file):
    dictbr={}
    file=open(mO_file,'r')
    f=file.readlines()
    file.close()
    for line in f:
        #REGEX ^:begin \d: digit, \s: space, \n: new line, .*: any set, (): group
        br=re.search(r'(^\d\.\d+[Ee][+\-]\d+)\s(.*)\n',line) #branchings
        if br:
            #groups: 2:(.*)-> decay name; 1: (^\d\.\d+[Ee][+\-]\d+)-> value  
            dictbr[br.group(2)]=float(br.group(1))

        tw=re.search(r'(^H\->2\*x).*(\d\.\d+[Ee][+\-]\d+).*\n',line) #totalwidth
        if tw:
            dictbr[tw.group(1)]=float(tw.group(2))

    return dictbr

if __name__ == '__main__':
    '''micrOMEGAS output in kk'''
    BR=read_micrOMEGAS_output('kk')
    totalwidth=BR['H->2*x']
    print totalwidth,BR['H -> s,S']
    
