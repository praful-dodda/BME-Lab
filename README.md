# BME-Lab
This repository is to share codes within BME labmates

The list of codes included and their functionalities:
<li> cmaq_camp_correction.m - 'For CAMP correction of CMAQ data'

2/8/23

Note: this forces from the middle decile. For example, if there are 10 bins, the function will make the slope increasing from the 5th bin forward, and decreasing from the 5th bin backwards.
  
<li> lambdaMonotonicCorr.m: exact function that Nora used in her CAMP correction to force curve monotonic. 
<li> monotonicCorrTest.m: generalized function for Nora's monotonic CAMP correction; does the exact same thing as lambdaMonotonicCorr but with better function syntax and broader variable names
