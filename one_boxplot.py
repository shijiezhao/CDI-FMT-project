"""
@Author: Shijie Zhao; @Sep/26/15; @Cambridge, MA
@Boxplot module
"""
import pylab
from zero_preprocessing import *
## Import for plotly 
import plotly.plotly as py
from plotly.graph_objs import *
py.sign_in('Python-Demo-Account', 'gwt101uhh0')
import numpy as np
import scipy.stats as ss

#######################################################################################################

## Construct df for different categories
Post = SubSetDataFrame(df=SDIs, label='type', value=['post'],meta_dict=meta_dict)
Pre = SubSetDataFrame(df=SDIs, label='type', value=['pre'],meta_dict=meta_dict)
Post_suc = SubSetDataFrame(df=Post, label='ssx_resolved', value=['Y'],meta_dict=meta_dict)
Post_fail = SubSetDataFrame(df=Post, label='ssx_resolved', value=['N'],meta_dict=meta_dict)
Donor = SubSetDataFrame(df=SDIs, label='type', value=['donor'],meta_dict=meta_dict)

## Save the files of p-value
y1 = Donor.values
y2 = Pre.values
y3 = Post_suc.values
y4 = Post_fail.values
y5 = Post.values

f = open('SDI_Pval.txt','w')
s = 'donor-pre: '
p = ss.mannwhitneyu(y1,y2)
S = s+str(p)+'\n'
f.write(S)
print s, p

s = 'donor-post: '
p = ss.mannwhitneyu(y1,y5)
S = s+str(p)+'\n'
f.write(S)
print s, p

s = 'pre-post: '
p = ss.mannwhitneyu(y2,y5)
S = s+str(p)+'\n'
f.write(S)
print s, p

s = 'suc-fail: '
p = ss.mannwhitneyu(y3,y4)
S = s+str(p)+'\n'
f.write(S)
print s, p


## Plot the boxplot
Zs = []
Zs.append(y1)
Zs.append(y2)
Zs.append(y3)
Zs.append(y4)
traces = [Box(y=z, boxpoints='all', pointpos=0) for z in Zs]
data = Data(traces)
py.plot(data, filename='basic-box-plot')