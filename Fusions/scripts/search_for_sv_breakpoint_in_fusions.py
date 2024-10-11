import sys, os
import argparse
import vcf

class Fusion: 
    def __init__(self, myline):
        data=myline.strip().split('\t')
        self.FusionName=data[0]
        self.JunctionReadCount=int(data[1])
        self.SpanningFragCount=int(data[2])
        self.est_J=float(data[3])
        self.est_S=float(data[4])
        self.SpliceType=data[5]
        self.LeftGene=data[6]
        self.LeftBreakpoint=data[7]
        self.RightGene=data[8]
        self.RightBreakpoint=data[9]
        self.JunctionReads=data[10]
        self.SpanningFrags=data[11]
        self.LargeAnchorSupport=data[12]
        self.FFPM=float(data[13])
        leftCoords=self.LeftBreakpoint.split(':')
        self.LeftChr=leftCoords[0]
        self.LeftLoc=int(leftCoords[1])
        rightCoords=self.RightBreakpoint.split(':')
        self.RightChr=rightCoords[0]
        self.RightLoc=int(rightCoords[1])
    def __str__(self):
        return "%s\n" % ("\t".join(data))
    def outstr(self):
        return "%s\t%d\t%d\t%s\t%s" % (self.FusionName, self.JunctionReadCount, self.SpanningFragCount, self.LeftBreakpoint, self.RightBreakpoint)

class FusionOutJunction:
    def __init__(self):
        self.Pos=0
        self.Qual=0
        self.Filter=""
        self.Type="."
        self.SubType="."
        
    def __str__(self):
        if self.Pos != 0 and len(self.Filter)==0:
            self.Filter="PASS"
        return ":".join(map(str, [self.Pos, self.Qual, self.Filter, self.Type]))
    
def get_best_SVbreakpoint(myfusOut, myChr, myLoc, vcf_reader):
    mydist=2*slop
    #sys.stderr.write("fetching: %s:%d-%d\n" % (myChr, myLoc-slop, myLoc+slop))
    recQUAL=0
    try:
        for rec in vcf_reader.fetch(myChr, myLoc-slop, myLoc+slop):
            if 'END' in rec.INFO.keys():
                recEND=rec.INFO['END']
            else: 
                recEND=rec.POS
            d=min(abs(rec.POS - myLoc), abs(recEND - myLoc))
            recQUAL=0
            if rec.QUAL is None:
                if ('SU' in rec.INFO.keys()):
                    recQAUL=rec.INFO['SU']
                elif ('SOMATICSCORE' in rec.INFO.keys()):
                    recQUAL=rec.INFO['SOMATICSCORE']
            else:
                recQUAL=rec.QUAL
            if (d < 2*slop and (recQUAL > myfusOut.Qual or d < mydist)):
                if (abs(rec.POS - myLoc) <= abs(recEND - myLoc)):
                    myfusOut.Pos=rec.POS
                else:
                    myfusOut.Pos=recEND
                myfusOut.Qual=recQUAL
                if ((rec.FILTER is None) or (len(rec.FILTER)==0)): 
                    if ('FT' in rec.samples[0].data._fields):
                        myfusOut.Filter=rec.samples[0].data.FT
                    else:
                        myfusOut.Filter="PASS"
                else:
                    myfusOut.Filter=",".join(rec.FILTER)
                myfusOut.SubType=rec.var_subtype
                myfusOut.Type=rec.INFO['SVTYPE']
                mydist=d
    except ValueError:
        pass

slop=1000

def find_svbreakpoints_in_fusions(fusionfile, svvcf, output): 
    outfh = open(output, 'w')
    vcf_reader=vcf.Reader(filename=svvcf)
    fusionfh = open(fusionfile, 'r')
    headline=fusionfh.readline() # Skip the first line
    for fline in fusionfh:
        myfusion=Fusion(fline)
        myfusLeft=FusionOutJunction()
        get_best_SVbreakpoint(myfusLeft, myfusion.LeftChr, myfusion.LeftLoc, vcf_reader)
        myfusRight=FusionOutJunction()
        get_best_SVbreakpoint(myfusRight, myfusion.RightChr, myfusion.RightLoc, vcf_reader)
        outline="%s\t%s\t%s\n" % (myfusion.outstr(), str(myfusLeft), str(myfusRight)) 
        outfh.write(outline)


def main(): 
    parser=argparse.ArgumentParser(description="merge SV breakpoints with gene fusion breakpoint")
    parser.add_argument('fusions', help='a star-fusion.fusion_candidates.preliminary* file')
    parser.add_argument('svvcf', help='a .vcf file with SV calls (either from MANTA, LUMPY, or GRIDSS, may add others')
    parser.add_argument('output', help='the name of the output file')
    args=parser.parse_args()
    find_svbreakpoints_in_fusions(args.fusions, args.svvcf, args.output)

if __name__ == "__main__":
    main()
