# Slepian_ocean_add-on

An add-on package to the [Slepian function Software](https://geoweb.princeton.edu/people/simons/software.html) for calculating the sea level anomaly in a regional ocean. Additional data procedures (e.g., recovering GAD and correcting Inverted barometer) were applied with a selective Gaussion smooth, proved effective in reducing the north-south stripings in GRACE products.

## Notice before use
1. This code is based on the Free Software from the Simons Laboratories (https://geoweb.princeton.edu/people/simons/software.html) with some changes for oceanic application (e.g., GAD, GIA, IB correction).
2. It is strongly recommended to familiarize yourself with the foundational work before using this add-on code.

**Required software:**<br>
[slepian_alpha](https://github.com/csdms-contrib/slepian_alpha)  
[slepian_bravo](https://github.com/csdms-contrib/slepian_bravo)  
[slepian_delta](https://github.com/csdms-contrib/slepian_delta)  
[slepian_zero](https://github.com/csdms-contrib/slepian_zero) (actually only guyotphysics.m)   
[m_map](https://www.eoas.ubc.ca/~rich/map.html)  

## Citation Information:
Please cite our work and the foundational work by <a href="https://polarice.geo.arizona.edu/">C. Harig</a> &amp; <a href="http://www.frederik.net">F. J. Simons</a> as appropriate if you find this package useful.  

Harig, Christopher and Frederik J. Simons. 
Mapping Greeenland's mass loss in space and time.
<i>Proc. Natl. Acad. Sc.</i>, 109(49), 19934-19937, 2012.
<a href="http://dx.doi.org/10.1073/pnas.1206785109">doi:10.1073/pnas.1206785109</a>

Zhongtian, Ma, Hok Sum Fok, Robert Tenzer, Jianli Chen. 
A Novel Slepian Approach for Determining Mass-term Sea Level from GRACE over the South China Sea. 

## Usage of this package
1. Download the required software (or codes) and this codes as mentioned above.
2. Run SCS_case.m for a case study in the South China Sea. The procedures of this case generally follow the guideline in [slepian_delta](https://github.com/csdms-contrib/slepian_delta).
