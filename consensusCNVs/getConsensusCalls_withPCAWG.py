import sys
import viola

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])

consensus = viola.read_vcf(sys.argv[1],variant_caller="manta")
#pcawg = viola.read_vcf(sys.argv[2],variant_caller='delly')

#merged_vcf = viola.merge([consensus, pcawg], threshold=100)

#out = [sys.argv[3],'.PCAWG.consensus.100bp.vcf']
#merged_vcf.to_vcf("".join(out))
