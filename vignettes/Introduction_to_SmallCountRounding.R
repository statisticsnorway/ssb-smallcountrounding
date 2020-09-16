## ----comment=NA, tidy = TRUE, include =FALSE----------------------------------
library(knitr) 
library(kableExtra)
library(SmallCountRounding)

cell_background = function(x,row,col,background){
  backGround = rep("white", 100)
  backGround[row] = background
  suppressWarnings(column_spec(x,col, background = backGround))
}
yellow = "#FFFF88"
green = "#F0FFE9"
green2 = "#88FF88"
z <- SmallCountData("exPSD")
set.seed(12345) 
a <- PLSrounding(z, "freq", 5)
k <- PLS2way(a, "original") 
ka <- PLS2way(a)
set.seed(12345) 
b <- PLSrounding(z, "freq", 5, formula = ~rows + cols)
kb <- PLS2way(b)


e6  <-  SmallCountData("e6")

eDimList <- SmallCountData("eDimList")

set.seed(12345)
e6a <-  PLSrounding(e6, "freq", 5)
set.seed(12345)
e6b <-  PLSrounding(e6, "freq", 5, formula = ~eu * year + geo * year)
set.seed(12345)
e6c <-  PLSrounding(e6[, -2], "freq", 5, hierarchies = eDimList)
set.seed(12345)
e6d <-  PLSrounding(e6[, -2], "freq", 5, hierarchies = eDimList, formula = ~geo * year)




options(knitr.kable.NA = '')


## ----comment=NA, tidy = TRUE, echo=FALSE--------------------------------------

kable(k, "html", caption = "**Table 1**: Original data in tabular form with row and column totals") %>%
  kable_styling(full_width = F, bootstrap_options = c("bordered"),font_size = 16,  position = "left")  %>%
  add_indent(1:4,1.6,all_cols =TRUE) %>%
  column_spec(1, bold = T,background = green)  %>%
  #cell_background(2,2,"#FFFF33") %>%
  #cell_background(2,5,"orange") %>%
  #column_spec(3, background = c("#FFFFFF","#FFFFFF","#FF1111","white"))  %>%
  column_spec(7, bold = T)  %>%
  row_spec(0, bold = T,background = green)  %>%
  row_spec(4, bold = T) 

## ----comment=NA, tidy = TRUE, echo=FALSE--------------------------------------

kable(ka, "html", caption = '**Table 2**: All small inner cell values (1-4) are rounded using 5 as rounding base.') %>%
  kable_styling(full_width = F, bootstrap_options = c("bordered"),font_size = 16,  position = "left")  %>%
  add_indent(1:4,1.5,all_cols =TRUE) %>%
  column_spec(1, bold = T,background = green)  %>%
  cell_background(2,2,yellow) %>%
  cell_background(2:3,3,yellow) %>%
  cell_background(1:3,4,yellow) %>%
  cell_background(1:2,5,yellow) %>%
  cell_background(1:3,6,yellow) %>%
  #cell_background(2,5,"orange") %>%
  #column_spec(3, background = c("#FFFFFF","#FFFFFF","#FF1111","white"))  %>%
  column_spec(7, bold = T)  %>%
  row_spec(0, bold = T,background = green)  %>%
  row_spec(4, bold = T) 

## ----comment=NA, tidy = TRUE, echo=FALSE--------------------------------------

kable(kb, "html", caption = '**Table 3**: Assuming only row and column totals to be published, necessary small inner cell values (1-4) are rounded using 5 as rounding base.') %>%
  kable_styling(full_width = F, bootstrap_options = c("bordered"),font_size = 16,  position = "left")  %>%
  add_indent(1:4,1.5,all_cols =TRUE) %>%
  column_spec(1, bold = T,background = green)  %>%
  cell_background(2:3,3,yellow) %>%
  cell_background(1:3,4,yellow) %>%
  cell_background(1:2,5,yellow) %>%
  cell_background(3,6,yellow) %>%
  #cell_background(2,5,"orange") %>%
  #column_spec(3, background = c("#FFFFFF","#FFFFFF","#FF1111","white"))  %>%
  column_spec(7, bold = T)  %>%
  row_spec(0, bold = T,background = green)  %>%
  row_spec(4, bold = T) 

## ----comment=NA, tidy = TRUE--------------------------------------------------
library(SmallCountRounding)
z <- SmallCountData("exPSD")
z

## ----eval=FALSE, tidy = TRUE--------------------------------------------------
#  a <- PLSrounding(z, freqVar = "freq", roundBase = 5)

## ----comment=NA, tidy = TRUE--------------------------------------------------
a$inner

## ----comment=NA, tidy = TRUE--------------------------------------------------
a$publish

## ----eval=FALSE, tidy = TRUE--------------------------------------------------
#  b <- PLSrounding(z, "freq", 5, formula = ~rows + cols)

## ----comment=NA, tidy = TRUE--------------------------------------------------
b$inner

## ----comment=NA, tidy = TRUE--------------------------------------------------
b$publish

## ----comment=NA, tidy = TRUE--------------------------------------------------
a
b

## ----comment=NA, tidy = FALSE-------------------------------------------------
f <- b$publish$original
g <- b$publish$rounded
print(c(
  maxdiff        =  max(abs(g - f)), 
  HDutility      =  HDutility(f, g), 
  meanAbsDiff    =  mean(abs(g - f)), 
  rootMeanSquare =  sqrt(mean((g - f)^2))
))

## ----comment=NA, tidy = FALSE-------------------------------------------------
summary(b)

## ----comment=NA, tidy = TRUE, echo=FALSE--------------------------------------
kable(e6, "html", caption = "**Table 4**: Input data") %>%
  kable_styling(full_width = F, bootstrap_options = c("bordered"),font_size = 14,  position = "left") %>%
  add_indent(1:6,0.4,all_cols =TRUE)

## ----comment=NA, tidy = TRUE, echo=FALSE--------------------------------------
kable(SmallCountData("eDimList")$geo, "html", caption = "**Table 5**:  Hierarchy, `geo`") %>%
  kable_styling(full_width = F, bootstrap_options = c("bordered"),font_size = 14,  position = "left")  %>%
  add_indent(1:6,1.2,all_cols =TRUE)

## ----comment=NA, tidy = TRUE, echo=FALSE--------------------------------------
kable(e6a$publish, "html", caption = "**Table 6**: Ouput data (publish)") %>%
  kable_styling(full_width = F, bootstrap_options = c("bordered"),font_size = 14,  position = "left") %>%
  column_spec(1:2, background = green) %>%
  cell_background( Match(e6a$inner,e6a$publish),4,c(yellow,green2)[1+(e6a$inner$difference==0)]) 
#cell_background( Match(e6a$inner[e6a$inner$difference!=0 , ],e6a$publish),4,yellow) %>%
#cell_background( Match(e6a$inner[e6a$inner$difference==0 , ],e6a$publish),4,green2) 

## ----comment=NA, tidy = TRUE, eval = TRUE-------------------------------------
e6  <-  SmallCountData("e6")             # As Table 4 
eDimList <- SmallCountData("eDimList")
eDimList

## ----comment=NA, tidy = FALSE, eval = FALSE-----------------------------------
#  PLSrounding(e6, "freq", 5)                                                      # a)
#  PLSrounding(e6, "freq", 5, formula = ~eu * year + geo * year)                   # b)
#  PLSrounding(e6[, -2], "freq", 5, hierarchies = eDimList)                        # c)
#  PLSrounding(e6[, -2], "freq", 5, hierarchies = eDimList, formula = ~geo * year) # d)

## ----comment=NA, tidy = FALSE, eval = TRUE------------------------------------
out <- PLSrounding(e6[-1, ], "freq", 5, removeEmpty = TRUE, inputInOutput = c(FALSE,TRUE))
out
out$inner
out$publish

