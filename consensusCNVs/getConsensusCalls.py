import sys
import viola

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])

manta = viola.read_vcf(sys.argv[1],variant_caller="manta")
gridss = viola.read_vcf(sys.argv[2],variant_caller='gridss')

merged_vcf = viola.merge([manta, gridss], threshold=150)
filtered_vcf = merged_vcf.filter('supportingcallercount >= 2')

out = [sys.argv[3],'.consensus.150bp.vcf']
filtered_vcf.to_vcf("".join(out))
