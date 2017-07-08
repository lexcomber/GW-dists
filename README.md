# GW distance metric paper

##Distance metric choice can both reduce and induce collinearity in geographically weighted regression

Alexis Comber (1), Binbin Lu (2) and Paul Harris (3)

1 School of Geography, University of Leeds, Leeds, LS2 9JT, UK \
2 School of Remote Sensing and Information Engineering, Wuhan University, Wuhan, China \
3 Rothamsted Research, North Wyke, EX20 2SB, UK |

# Abstract
This paper explores the impact of distance metric on collinearity on local regression models such as Geographically Weighted Regression (GWR). Using a case study of data related to house price in Hanoi Vietnam, and by fully varying both power and rotation to create different Minkowski distances, the analysis shows that local collinearity can be both negatively affected by distance metric choice. The Minkowski distance that maximised collinearity was approximate to Manhattan distance with (power = 0.70) with a rotation of 30 degrees and that which minimised collinearity was specified with power = 0.05 and a rotation of 70 degrees. The discussion considers the exploration and selection of different distance metrics to minimise collinearity as an as alternative to elastic nets, lassos and ridge regression in local regression approaches. It suggests potential areas of further work by extending the approach to include measures of model fit in the manner of Lu et al. (2016) Use these methods to select distance metric to minimise collinearity (i.e. maximizing the proportion of the proportion of CNs < 30 could be another criterion for distance metric choice) and the need to investigate the relationship between kernel bandwidth, distance metrics and collinearity. 

