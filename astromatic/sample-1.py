
import sextractor

# Create a SExtractor instance
sex = sextractor.SExtractor()

# Modify the SExtractor configuration
sex.config['GAIN'] = 0.938
sex.config['PIXEL_SCALE'] = .19
sex.config['VERBOSE_TYPE'] = "FULL"
sex.config['CHECKIMAGE_TYPE'] = "BACKGROUND"

# Add a parameter to the parameter list
sex.config['PARAMETERS_LIST'].append('FLUX_BEST')

# Lauch SExtractor on a FITS file
sex.run("lm00310SbMcl2426.fits")

# Read the resulting catalog [first method, whole catalog at once]
catalog = sex.catalog()
for star in catalog:
    print star['FLUX_BEST'], star['FLAGS']
    if (star['FLAGS'] & sextractor.BLENDED):
        print "This star is BLENDED"

# Read the resulting catalog [second method, whole catalog at once]
catalog_name = sex.config['CATALOG_NAME']
catalog_f = sextractor.open(catalog_name)
catalog = catalog_f.readlines()
for star in catalog:
    print star['FLUX_BEST'], star['FLAGS']
    if (star['FLAGS'] & sextractor.BLENDED):
        print "This star is BLENDED"

catalog_f.close()
        
# Read the resulting catalog [third method, star by star]
catalog_name = sex.config['CATALOG_NAME']
catalog_f = sextractor.open(catalog_name)
star = catalog_f.readline()
while star:
    print star['FLUX_BEST'], star['FLAGS']
    if (star['FLAGS'] & sextractor.BLENDED):
        print "This star is BLENDED"
    star = catalog_f.readline()

catalog_f.close()

# Removing the configuration files, the catalog and
# the check image
sex.clean(config=True, catalog=True, check=True)

