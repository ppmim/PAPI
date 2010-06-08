
# Through sexcatalog module
import sexcatalog

# Read a SExtractor ASCII catalog

catalog_name = "lm038.cat"

# First method: read the whole catalog at once
catalog_f = sexcatalog.open(catalog_name)
catalog = catalog_f.readlines()
for star in catalog:
    print star['FLUX_BEST'], star['FLAGS']
    if (star['FLAGS'] & sexcatalog.BLENDED):
        print "This star is BLENDED"
catalog_f.close()

# Second method: read the catalog star by star
catalog_f = sexcatalog.open(catalog_name)
for star in catalog_f:
    print star['FLUX_BEST'], star['FLAGS']
    if (star['FLAGS'] & sexcatalog.BLENDED):
        print "This star is BLENDED"
catalog_f.close()

