# Slepian_ocean_add-on

An add-on package to the [Slepian function Software](https://geoweb.princeton.edu/people/simons/software.html) for calculating the sea level variations in a regional ocean. Additional data procedures (e.g., recovering GAD and correcting Inverted barometer) are applied with a selective Gaussion smooth, which proves effective in reducing the north-south stripings in GRACE products.

## Notice before use
1. This code is based on the [foundational Software from the Simons Laboratories](https://geoweb.princeton.edu/people/simons/software.html) with some changes for ocean application (e.g., GAD, GIA, IB correction).
2. It is strongly recommended to familiarize yourself with the foundational work by <a href="https://polarice.geo.arizona.edu/">C. Harig</a> &amp; <a href="http://www.frederik.net">F. J. Simons</a> before using this add-on code.

**Required software:**<br>
[slepian_alpha](https://github.com/csdms-contrib/slepian_alpha)  
[slepian_bravo](https://github.com/csdms-contrib/slepian_bravo)  
[slepian_delta](https://github.com/csdms-contrib/slepian_delta)  
[slepian_zero](https://github.com/csdms-contrib/slepian_zero) (actually only guyotphysics.m)   
[m_map](https://www.eoas.ubc.ca/~rich/map.html)  

## Citation Information:
Please cite our work and the foundational work by C. Harig and F. J. Simons as appropriate if you find this package useful.  

Harig, Christopher and Frederik J. Simons. 
Mapping Greenland's mass loss in space and time.
<i>Proc. Natl. Acad. Sc.</i>, 109(49), 19934-19937, 2012.
doi:https://doi.org/10.1073/pnas.1206785109

Zhongtian, Ma, Hok Sum Fok, Robert Tenzer, Jianli Chen. 
A Novel Slepian Approach for Determining Mass-term Sea Level from GRACE over the South China Sea. <i>Int. J. Appl. Earth. Obs. Geoinf.</i>, 132, 104065, 2024.
doi:https://doi.org/10.1016/j.jag.2024.104065

## Usage of this package
1. Download the required software (or codes) and this codes as mentioned above.
2. Download the Spherical Harmonic coefficients and the associated correction files (e.g., degree 1, degree 20 correction) from official websites.
3. Run SCS_case.m for a case study in the South China Sea. The procedures of this case generally follow the guideline in [slepian_delta](https://github.com/csdms-contrib/slepian_delta).
4. (Optional) Compare this mass-term sea level from Slepian functions with other datasets (e.g., Spherical harmonics, Mascon) in the South China Sea. You can download the mass-term sea level from other datasets in [Zenodo](https://zenodo.org/records/12684255).  
