# GW distance metric paper - code and data
This analysis is described in a paper soon to be submitted. The code and data are all described in the file `GWDists_paper_code.R`. Some of the file paths (`setwd()`) may have to be changed to reflect your local directory.

## Paper Title: Distance metric choice can both reduce and induce collinearity in geographically weighted regression

Alexis Comber<sup>1</sup>, Khanh Chi<sup>2</sup>, Man Quang Huy<sup>3</sup>, Quan Nguyen<sup>4</sup>, Binbin Lu<sup>5</sup>, Hoang Huu Phe<sup>6</sup> and Paul Harris<sup>7</sup> 

<sup>1</sup>School of Geography, University of Leeds, LS2 9JT, UK\
<sup>2</sup>GeoViet Consulting Co., Ltd, Hanoi, Vietnam\
<sup>3</sup>Vietnam National University Hanoi\
<sup>4</sup>National University of Civil Engineering, Hanoi, Vietnam\
<sup>5</sup>Wuhan University, Wuhan, China\
<sup>6</sup>Vinaconex R&D, Hanoi, Vietnam\
<sup>7</sup>Rothamsted Research, North Wyke, EX20 2SB, UK\

# Abstract
This paper explores the impact of different distance metrics on collinearity in local regression models such as Geographically Weighted Regression (GWR). Using a case study of house price data collected in Hà Nội, Vietnam, and by fully varying both power and rotation parameters to create different Minkowski distances, the analysis shows that local collinearity can be both negatively and positively affected by distance metric choice. The Minkowski distance that maximised collinearity in GWR was approximate to a Manhattan distance with (power = 0.70) with a rotation of 30, and that which minimised collinearity was parameterised with power = 0.05 and a rotation of 70. The results indicate that distance metric choice can provide a useful extra tuning component to address local collinearity issues in spatially varying coefficient modelling and that understanding the interaction of distance metric and collinearity can provide insight into the nature and structure of the data relationships. The discussion considers first, the exploration and selection of different distance metrics to minimise collinearity as an alternative to localised ridge regression, lasso and elastic net approaches. Second, it discusses the how distance metric choice could extend the methods that additionally optimise local model fit (lasso and elastic net) by selecting a distance metric that further helped minimise local collinearity. Third, it identifies the need to investigate the relationship between kernel bandwidth, distance metrics and collinearity as an area of further work.
**Keywords**: GWR; distance metrics; model fit; collinearity.

November 2017
