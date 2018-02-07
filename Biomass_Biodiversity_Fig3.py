## Code adapted from https://github.com/LennonLab/ScalingMicroBiodiversity/Fig1.py by CM

from __future__ import division
import  matplotlib.pyplot as plt

import pandas as pd
import linecache
import numpy as np
import os
import sys
import scipy.stats as stats

import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table


mydir = os.getcwd()


## FIGURES A and B Make read-based N subsurface dataframe and fit
def LL_plot_ci_pi(x,y,fillcolor,edgecolor,transparency,includePI,includeCI):
	# Build dataframe from x,y
	dt = pd.DataFrame({'N': x}) 
	dt['y'] = y

	# make fit, recover r-squared, a, and b and set up summary table
	regfit = smf.ols('y ~ N',dt).fit() 
	r2 = regfit.rsquared
	params = regfit.params
	a = params[0]
	b = params[1]
	st, data, ss2 = summary_table(regfit, alpha =0.05)

	#pull confidence and prediction intervals
	fittedvalues = data[:,2]
	predict_mean_se  = data[:,3]
	predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
	predict_ci_low, predict_ci_upp = data[:,6:8].T	

	# build data_intervals object to store 95% PIs and CIs
	data_intervals = {'N': x, 'predict_low': predict_ci_low, 'predict_upp': predict_ci_upp, 'conf_low': predict_mean_ci_low, 'conf_high': predict_mean_ci_upp}
	df_intervals = pd.DataFrame(data=data_intervals)
	df_intervals_sort = df_intervals.sort_values(by='N')
		
	#Make scatter plot and best fit line
	plt.scatter(x,y,color = fillcolor, alpha= transparency , s = 8, linewidths=0.5, edgecolor=edgecolor)
	plt.plot(x, fittedvalues, '-', lw=2, color = fillcolor)
	
	#Draw PIs and CIs
	if 'yes' in includePI:
		draw_pred_interval(df_intervals_sort['N'], df_intervals_sort['predict_low'], df_intervals_sort['predict_upp'],fillcolor)				
	if 'yes' in includeCI:
		plt.fill_between(df_intervals_sort['N'], df_intervals_sort['conf_high'], df_intervals_sort['conf_low'], color=fillcolor, alpha = 0.5,lw=0.0)

	return r2, a, b

## FIGURE C species vs volume
def sar_ci_pi(dt,type, fillcolor,edgecolor,transparency,includePI,includeCI,md,legendtitle):
	#Build dataframe from imported dt object
	dt['x']= np.log10(dt['volume.filtered.in.L'])
	dt['y'] = np.log10(dt[type])
	
	# make fit, recover r-squared, a, and b and set up summary table
	regfit = smf.ols('y ~ x',dt).fit()
	r2 = regfit.rsquared
	params = regfit.params
	a=params[0]
	b = params[1]
	st, data, ss2 = summary_table(regfit, alpha =0.05)

	#pull confidence and prediction intervals
	fittedvalues = data[:,2]
	predict_mean_se  = data[:,3]
	predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
	predict_ci_low, predict_ci_upp = data[:,6:8].T	

	# build data_intervals object to store 95% PIs and CIs
	data_intervals = {'x': dt['x'], 'predict_low': predict_ci_low, 'predict_upp': predict_ci_upp, 'conf_low': predict_mean_ci_low, 'conf_high': predict_mean_ci_upp}
	df_intervals = pd.DataFrame(data=data_intervals)
	df_intervals_sort = df_intervals.sort_values(by='x')
	
	# add points based on primer
	for kind in md:
		d = dt.loc[dt['primer']==kind]
		plt.scatter(d.x, d.y,marker = md[kind],color = fillcolor, alpha= transparency , s = 25, linewidths=0.5, edgecolor=edgecolor,label=None)		

	# add best fit line
	plt.plot(dt['x'], fittedvalues, '-', lw=2, color = fillcolor, label=legendtitle)
	
	#Draw PIs and CIs
	if 'yes' in includePI:
		draw_pred_interval(df_intervals_sort['x'], df_intervals_sort['predict_low'], df_intervals_sort['predict_upp'],fillcolor)

	if 'yes' in includeCI:
		plt.fill_between(df_intervals_sort['x'], df_intervals_sort['conf_high'], df_intervals_sort['conf_low'], color=fillcolor, alpha = 0.5,lw=0.0)

	return r2, a, b


def dd_plot_ci_pi(dt,fillcolor,edgecolor,transparency,includeCI,md):
# 	Construct data frame
	dt['D'] = np.log10(dt['x'])
	dt['y'] = dt['similarity']
	
	# make fit, recover r-squared, a, and b and set up summary table
	regfit = smf.ols('y ~ D',dt).fit()
	r2 = regfit.rsquared
	params = regfit.params
	a=params[0]
	b = params[1]
	st, data, ss2 = summary_table(regfit, alpha =0.05)

	#pull confidence and prediction intervals
	fittedvalues = data[:,2]
	predict_mean_se  = data[:,3]
	predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
	predict_ci_low, predict_ci_upp = data[:,6:8].T	

	# build data_intervals object to store 95% PIs and CIs
	data_intervals = {'D':  dt['D'], 'predict_low': predict_ci_low, 'predict_upp': predict_ci_upp, 'conf_low': predict_mean_ci_low, 'conf_high': predict_mean_ci_upp}
	df_intervals = pd.DataFrame(data=data_intervals)
	df_intervals_sort = df_intervals.sort_values(by='D')

	# Make scatter plots
	for kind in md:
		print kind
		d = dt.loc[dt['primer']==kind]
		if kind == 'v6':
			plt.scatter(d.D, d.y,marker = md[kind],color = 'grey',edgecolor = 'crimson', alpha= transparency , s = 25, linewidths=0.5,label=None)
 		else:
 			plt.scatter(d.D, d.y,marker = md[kind],color = 'maroon',edgecolor = 'crimson', alpha= transparency , s = 25, linewidths=0.5,label=None)		
		
	plt.plot(dt['D'], fittedvalues, '-', lw=2, color = fillcolor,label=r'All; $R^2=0.37$')
	
	if 'yes' in includeCI:
		plt.fill_between(df_intervals_sort['D'], df_intervals_sort['conf_high'], df_intervals_sort['conf_low'], color=fillcolor, alpha = 0.5,lw=0.0)
	print r2
	return r2, a, b

def special_plot_ci_pi(dt,variabletype,var,color,linetype,transparency,includeCI,md):
# 	Construct data frame
	dt['N'] = np.log10(dt['x'])
	dt['y'] = dt['similarity']
	
	# pull site of interest
	dt = dt.loc[dt[variabletype]==var]
	
	# make fit, recover r-squared, a, and b and set up summary table
	regfit = smf.ols('y ~ N',dt).fit()
	r2 = regfit.rsquared
	params = regfit.params
	a=params[0]
	b = params[1]
	st, data, ss2 = summary_table(regfit, alpha =0.05)

	#pull confidence and prediction intervals
	fittedvalues = data[:,2]
	predict_mean_se  = data[:,3]
	predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
	predict_ci_low, predict_ci_upp = data[:,6:8].T	

	# build data_intervals object to store 95% PIs and CIs
	data_intervals = {'N':  dt['N'], 'predict_low': predict_ci_low, 'predict_upp': predict_ci_upp, 'conf_low': predict_mean_ci_low, 'conf_high': predict_mean_ci_upp}
	df_intervals = pd.DataFrame(data=data_intervals)
	df_intervals_sort = df_intervals.sort_values(by='N')
	
	# plot 
	legendtitle = r'all data; $R^2=0.35$'.replace('all data',var).replace('0.35',str(round(r2,2)))
	if var == 'v6':
		legendtitle = r'all data; $R^2=0.35$'.replace('all data',var).replace('0.35','0.00')

	plt.plot(dt['N'], fittedvalues, linetype, lw=1, color = color,label=legendtitle)

	if var == 'USA':
		for kind in md:
			print kind
			d = dt.loc[dt['primer']==kind]
			plt.scatter(d.N, d.y,marker = md[kind],color = 'white',edgecolor = 'white', alpha= 1 , s = 25, linewidths=0.5,label=None)				
			plt.scatter(d.N, d.y,marker = md[kind],color = color,edgecolor = 'crimson', alpha= transparency , s = 25, linewidths=0.5,label=None)	

	if 'yes' in includeCI:
		plt.fill_between(df_intervals_sort['N'], df_intervals_sort['conf_high'], df_intervals_sort['conf_low'], color=color, alpha = 0.25,lw=0.0)

	print r2
	return r2, a, b

def draw_pred_interval(x,lower,upper, fillcolor):
	plt.plot(x, lower, '--', lw=1, color = fillcolor)
	plt.plot(x, upper, '--', lw=1, color = fillcolor)

def Fig3(condition, ones):

## Declare whether or not to use singleton data
    tail = str()
    if ones is False:
        tail = '-SADMetricData_NoMicrobe1s.txt'
    elif ones is True:
        tail = '-SADMetricData.txt'

    datasets = []
    GoodNames = []
    emp = str()

## Declare whether to use open or closed reference OTU clustering
    if condition == 'open':
    	emp = 'EMPopen'
    	subsurface = 'SUBopen'
    elif condition == 'closed': 
    	emp = 'EMPclosed'
    	subsurface = 'SUBclosed'

## GoodNames is a variable designating subdirectories to use
	GoodNames = [emp, 'HMP', 'BIGN', 'TARA', 'BOVINE', 'HUMAN', 'LAUB', 'SED', 'CHU', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO',  'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA']
	
## Add subsurface data N=reads
	GoodNames.append(subsurface)

## Add subsurface data N=cells	
	GoodNames.append('nSUBopen')
	GoodNames.append('nSUBclosed')
	
    for name in os.listdir(mydir +'/data/micro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'/data/micro/'+name+'/'+name+tail
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines


## Replace searching the macro folder with using subsurface data

    for name in os.listdir(mydir +'/data/subsurface'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'/data/subsurface/'+name+'/'+name+tail
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'subsurface', num_lines])
        print name, num_lines

    for name in os.listdir(mydir +'/data/nsubsurface'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'/data/nsubsurface/'+name+'/'+name+tail
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'nsubsurface', num_lines])
        print name, num_lines

    metrics = ['Rarity, '+r'$log_{10}$',
            'Dominance, '+r'$log_{10}$',
            'Evenness, ' +r'$log_{10}$',
            'Richness, ' +r'$log_{10}$',] #+r'$(S)^{2}$']

    fig = plt.figure(figsize=(11.69,8.27))
    fs = 8 # font size used across figures

    subplotdictionary = {1:0,3:1}
   
## Go through each index and generate a plot    
    for index, i in enumerate(metrics):
        if index in set([1,3]):        
            fig.add_subplot(2, 2, subplotdictionary[index]+1)
        metric = i

## Variables with Mac* now correspond to subsurface data, N=reads, being pulled in from the subsurface directory
        MicIntList, MicCoefList, MacIntList, MacCoefList, R2List, metlist = [[], [], [], [], [], []]


        Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

### ORIGINAL CODE used an iterative approach to obtain fits for micro data. It would select a subset of each micro subdirectory samples, assigning them to a variable numlines


        Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

        numMac = 0
        numMic = 0
        radDATA = []

## Build data structure to work with in matplotlib
        for dataset in datasets:

            name, kind, numlines = dataset

            lines = range(1, numlines+1) 
            
            path = mydir+'/data/'+kind+'/'+name+'/'+name+tail

            for line in lines:
                data = linecache.getline(path, line)
                radDATA.append(data)



## CM adjustment. Highlight mgp68 aka LAUB, HMP, subsurface, nsubsurface datasets

        HMPind = []
        LAUBind = []
        subsurfaceInd = []
        nsubsurfaceClosedInd = []


        for ti,data in enumerate(radDATA):

            data = data.split()
            name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

            N = float(N)
            S = float(S)

            Nlist.append(float(np.log10(N)))
            Slist.append(float(np.log10(S)))

            ESimplist.append(float(np.log10(float(ESimp))))

            KindList.append(kind)

            BPlist.append(float(BP))
            NmaxList.append(float(np.log10(float(Nmax))))

            # log-modulo transformation of skewnness
            lms = np.log10(np.abs(float(skew)) + 1)
            if skew < 0: lms = lms * -1
            rareSkews.append(float(lms))


            if name == 'HMP':
                HMPind.append(ti)
            if name == 'LAUB':
                LAUBind.append(ti)

            if kind == 'subsurface':
                subsurfaceInd.append(ti)
            if kind == 'nsubsurfaceclosed':
                nsubsurfaceClosedInd.append(ti)


        if index == 0: metlist = list(rareSkews)
        elif index == 1: metlist = list(NmaxList)
        elif index == 2: metlist = list(ESimplist)
        elif index == 3: metlist = list(Slist)

        
        ## Build dataframe
        d = pd.DataFrame({'N': list(Nlist)})
        d['y'] = list(metlist)
        d['Kind'] = list(KindList)

        ## PREPARE TO PLOT
        print index
        
        if index == 1:
            plt.ylim(0, 8)
            plt.xlim(0, 8.5)

        elif index == 3:
            plt.ylim(0.9, 5.0)
            plt.xlim(0, 12.5)


#####
##        ## Get 95% prediction interval of micro data
######

        ## Define micro subset
        MicListX = []
        MicListY = []

        for j, k in enumerate(KindList):
            if k == 'micro':
                MicListX.append(Nlist[j])
                MicListY.append(metlist[j])
                
        ## Get subsurface subset
        SUBd = d.loc[d['Kind'] == 'subsurface']
        SUBdX = [Nlist[i] for i in subsurfaceInd]
        SUBdY =[metlist[i] for i in subsurfaceInd]

        ## Get HMP data
        HMPListX = [Nlist[i] for i in HMPind]
        HMPListY = [metlist[i] for i in HMPind]
        
        ## Get mgp68 soil data
        LAUBListX = [Nlist[i] for i in LAUBind]
        LAUBListY = [metlist[i] for i in LAUBind]

## Make cell-based N subsurface dataframe and fit

        cClosedSUBd = d.loc[d['Kind'] == 'nsubsurfaceclosed']
        cClosedSUBdX = [Nlist[i] for i in nsubsurfaceClosedInd]
        cClosedSUBdY =[metlist[i] for i in nsubsurfaceClosedInd]

        if index == 3: # include Ncells for S~N
            cClosedSUBdR2, cClosedSUBdA, cClosedSUBdB = LL_plot_ci_pi(cClosedSUBdX,cClosedSUBdY,'maroon','Crimson',0.75,'yes','yes')

# plot everything 
        if index in set([1,3]):       
            plt.scatter(MicListX, MicListY, color = 'Grey', alpha= 0.15 , s = 8, linewidths=0.5, edgecolor='Steelblue')
            SUBdR2, SUBdA, SUBdB = LL_plot_ci_pi(SUBdX,SUBdY,'red','Crimson',0.8,'yes','yes')
            HMPR2, HMPA, HMPB = LL_plot_ci_pi(HMPListX, HMPListY,'darkblue','SteelBlue',0.6,'no','no')
            LAUBR2, LAUBA, LAUBB = LL_plot_ci_pi(LAUBListX, LAUBListY,'limegreen','SteelBlue',0.6,'no','no')
            print 'subsurface reads = %f*N^%f' %(round(10**SUBdA,2), round(SUBdB,2))
            print 'Soil = %f*N^%f' %(round(10**LAUBA,2), round(LAUBB,2))
            print 'HMP = %f*N^%f' %(round(10**HMPA,2), round(HMPB,2))
            plt.xlabel('$log$'+r'$_{10}$'+'($N$)', fontsize=fs+2)
            plt.ylabel(metric, fontsize=fs+2)
            plt.tick_params(axis='both', which='major', labelsize=fs-3)

# Plot fits from LL2016
        def graph(formula, x_range):  
            x = np.array(x_range)  
            y = eval(formula)
            plt.plot(x, y, ls = '-', lw=3, c='SteelBlue')  

        def graphS(formula, x_range):  
            x = np.array(x_range)
            y = eval(formula)
            plt.plot(x[:len(x_range)-1]  , y[1:], ls = '-', lw=3, c='SteelBlue')  


 # Make Fig A and B      
        if index == 0:
            print 'do nothing'

        elif index == 1:
            graph('np.log10(0.47*(10**x)**0.92)',range(13))
            if len(str(round(SUBdR2,2))) ==3:
                SUBdR2 = str(round(SUBdR2,2))+'0'
            plt.text(0.35, 7, r'$Subsurface N_{reads}$'+ ' $R^2$'+ '=' + SUBdR2, fontsize=fs, color='red')
            plt.text(0.35, 6,  r'$Soil$'+ ' $R^2$' + '=' + str(round(LAUBR2,2)), fontsize=fs, color='limegreen')
            plt.text(0.35, 5,  r'$HMP$'+ ' $R^2$' + '=' + str(round(HMPR2,2)), fontsize=fs, color='darkblue')

        elif index == 2:
            print 'do nothing'
            
        elif index == 3:
            graphS('np.log10(1.72*(10**x)**0.39)',range(13))

            print 'subsurface cells = %f*N^%f' %(round(10**cClosedSUBdA,2), round(cClosedSUBdB,2))

            plt.text(0.35, 4.5, r'$Subsurface N_{reads}$'+ ' $R^2$'+ '=' + str(round(SUBdR2,2)), fontsize=fs, color='red')
            plt.text(0.35, 4.0, r'$Subsurface N_{cells}$'+ ' $R^2$' + '=' + str(round(cClosedSUBdR2,2)), fontsize=fs, color='maroon')
            plt.text(0.35, 3.5,  r'$Soil$'+ ' $R^2$' + '=' + str(round(LAUBR2,2)), fontsize=fs, color='limegreen')
            plt.text(0.35, 3.0,  r'$HMP$'+ ' $R^2$' + '=' + str(round(HMPR2,2)), fontsize=fs, color='darkblue')




######## BUILD THE S'~Volume plot

    fig.add_subplot(2, 2, 3)
    plt.ylim(1, 5.5)
    plt.xlim(-1, 5.5)
    
    df = pd.read_csv(mydir + '/data/SubsurfaceVolumes/OTU_S_Ncells.csv')
    
    mkr_dict = {'v1v3': 'h', 'v3v5': '+', 'v4v5': '^', 'v6':'o','v6v4':'D'}
    ClosedSARR2, ClosedSARA, ClosedSARB= sar_ci_pi(df,'closedref','maroon','Crimson',0.75,'no','yes',mkr_dict,r'closed reference OTUs; $R^2=0.21$')
    ClosedSARR2, ClosedSARA, ClosedSARB= sar_ci_pi(df,'openref','lightcoral','Crimson',0.75,'no','yes',mkr_dict,r'open reference OTUs; $R^2=0.20$')

    plt.xlabel('Volume in L, $log$'+r'$_{10}$'+'($V$)', fontsize=fs+2)
    plt.ylabel('$log$'+r'$_{10}$'+"($S'$)", fontsize=fs+2)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.legend(fontsize=fs,loc='lower right')

######## BUILD THE distance decay plot

    fig.add_subplot(2, 2, 4)
    plt.ylim(0, 1)
    plt.xlim(-3.2, 5)
    
    df = pd.read_csv(mydir + '/data/SubsurfaceDistanceDecay/closedref_distance_sorensen_Depth_no_0.csv')
    
    
    mkr_dict = {'v1v3': 'h', 'v3v5': '+', 'v6':'o','v6v4':'D','different primers': '*', 'v4v5': '^'}
    colordict = {'FennoScandian Shield': 'cyan', 'different regions': 'black', 'USA': 'magenta', 'South Africa':'darkorange','Iceland':'darkgoldenrod'}

    markerList = []
    colorList = []

    for item in df['primer']:
	    markerList.append(mkr_dict[item])
    for item in df['locality']:
	    colorList.append(colordict[item])

    df['marker']=markerList
    df['color']=colorList

    dd_plot_ci_pi(df,'maroon','red',0.75,'yes',mkr_dict)
    special_plot_ci_pi(df,'primer','v6','Grey','-',0.9,'yes',mkr_dict)
    special_plot_ci_pi(df,'locality','USA','magenta','-',0.9,'yes',mkr_dict)
    
    
    plt.xlabel('Distance in km, $log$'+r'$_{10}$'+'(km)', fontsize=fs+2)
    plt.ylabel('Sorensen similarity', fontsize=fs+2)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.legend(fontsize=fs,loc='upper center',ncol=3)

    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    plt.show()
    plt.close()

    return



EMPcondition = ['closed']
Singletons = [False]

for condition in EMPcondition:
    for ones in Singletons:
        Fig3(condition, ones)
