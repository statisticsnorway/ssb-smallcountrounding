# Small Count Rounding of Tabular Data

A statistical disclosure control tool to protect frequency tables in
cases where small values are sensitive. The main function,
[`PLSrounding`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLSrounding.md),
performs small count rounding of necessary inner cells (Heldal, 2017) so
that all small frequencies of cross-classifications to be published
(publishable cells) are rounded. This is equivalent to changing micro
data since frequencies of unique combinations are changed. Thus,
additivity and consistency are guaranteed. This is performed by an
algorithm inspired by partial least squares regression (Langsrud and
Heldal, 2018).  

## References

Heldal, J. (2017): “The European Census Hub 2011 Hypercubes - Norwegian
SDC Experiences”. In: *Work Session on Statistical Data
Confidentiality*, Skopje, The former Yugoslav Republic of Macedonia,
September 20-22 , 2017.

Langsrud, Ø. and Heldal, J. (2018): “An Algorithm for Small Count
Rounding of Tabular Data”. Presented at: *Privacy in statistical
databases*, Valencia, Spain. September 26-28, 2018.
<https://www.researchgate.net/publication/327768398_An_Algorithm_for_Small_Count_Rounding_of_Tabular_Data>

## See also

Useful links:

- <https://github.com/statisticsnorway/ssb-smallcountrounding>

- <https://statisticsnorway.github.io/ssb-smallcountrounding/>

- Report bugs at
  <https://github.com/statisticsnorway/ssb-smallcountrounding/issues>

## Author

**Maintainer**: Øyvind Langsrud <oyl@ssb.no>

Authors:

- Øyvind Langsrud <oyl@ssb.no>

- Johan Heldal
