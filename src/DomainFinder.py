import argparse
import pyBigWig
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

parser = argparse.ArgumentParser()
parser.add_argument('--file_1',type=str, required=True, help='name of the first input file')
parser.add_argument('--file_2',type=str,required=True, help='name of the second input file')
parser.add_argument('--list_of_chrom',nargs='+', type=str, required=True, help='List of chromosome names separated by comma')
parser.add_argument('--bp_resolution_list', nargs='+', type=int, required=True, help='List of bp resolutions separated by comma')
args = parser.parse_args()
file_1 = args.file_1
file_2 = args.file_2
list_of_chrom = args.list_of_chrom[0].split(',')
bp_resolution_list = args.bp_resolution_list

print(file_1)
print(file_2)
print(list_of_chrom)
print(bp_resolution_list)

def find_consecutive_indexes(arr):
    consecutive_indexes = []

    if(len(arr) == 1):
        consecutive_indexes.append((0, 0))
        return consecutive_indexes
    flag = 0
    i_1 = 0
    i = 0
    for i in range(len(arr) - 1):
        if arr[i] + 1 == arr[i + 1]:
            if flag == 0:
                flag = 1
                i_1 = i
        if arr[i] + 1 != arr[i + 1]:
            if flag == 1:
                consecutive_indexes.append((i_1, i))
                i_1 = i+1
            if flag == 0:
                consecutive_indexes.append((i, i))
                i_1 = i+1
            flag = 0
           
    if flag == 1:
        consecutive_indexes.append((i_1, i+1))
    if flag == 0:
        consecutive_indexes.append((i+1, i+1))
    return consecutive_indexes
    
def find_domains(signal_list, bp_resolution):
    domains_list = []
    for signal in signal_list:
        sample_mean = np.mean(signal)
        sample_std = np.std(signal)
        neutral = sample_std

        sample_mean_plus = np.mean(np.take(signal,np.where(signal>=np.float64(0))))
        sample_std_plus = np.std(np.take(signal,np.where(signal>=np.float64(0))))

        sample_mean_minus = np.mean(np.take(signal,np.where(signal<np.float64(0))))
        sample_std_minus = np.std(np.take(signal,np.where(signal<np.float64(0))))

        ii_below = np.where(signal<sample_mean_minus-sample_std_minus)
        ii_above = np.where(signal>sample_mean_plus+sample_std_plus)
        ii_neutral = np.where((signal<=sample_mean_plus+sample_std_plus) & (signal>=sample_mean_minus-sample_std_minus))

        c = ['start', 'end']
        
        a = pd.DataFrame(find_consecutive_indexes(ii_below[0]),columns=c)
        if not a.empty and len(ii_below[0]) > 0:
            below_domains = pd.DataFrame({'start':ii_below[0][a.iloc[:, 0]],'end':ii_below[0][a.iloc[:, 1]],'type':"below",'color':"palevioletred"})
        else:
            below_domains = pd.DataFrame()

        a = pd.DataFrame(find_consecutive_indexes(ii_above[0]),columns=c)
        if not a.empty and len(ii_above[0]) > 0:
            above_domains = pd.DataFrame({'start':ii_above[0][a.iloc[:, 0]],'end':ii_above[0][a.iloc[:, 1]],'type':"above",'color':"cornflowerblue"})
        else:
            above_domains = pd.DataFrame()

        a = pd.DataFrame(find_consecutive_indexes(ii_neutral[0]),columns=c)
        if not a.empty and len(ii_neutral[0]) > 0:
            neutral_domains = pd.DataFrame({'start':ii_neutral[0][a.iloc[:, 0]],'end':ii_neutral[0][a.iloc[:, 1]],'type':"neutral",'color':"lightgrey"})
        else:
            neutral_domains = pd.DataFrame()

        domains = pd.concat([below_domains,above_domains,neutral_domains],ignore_index=True).sort_values(by='start',ascending=False)
        ##Get total signal in each domain
        domain_mean_signal = []
        for index, row in domains.iterrows():
            dd = round(sum(signal[row['start']:row['end']+1])/(row['end']-row['start']+1))
            domain_mean_signal.append(dd)
 
        domains["domain_mean_signal"] = domain_mean_signal
        domains_list.append(domains)
        
    return domains_list

def plot_domains_multiple_resolution(domains_list, chrom_names, chrom_sizes,bp_resolution_list):
    domains_list_with_names = zip(domains_list,chrom_names,chrom_sizes)
    for dd in domains_list_with_names:
        tick_step, x_label = get_x_tick_step(dd[2])
        tick_labels = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
        domain_bar_heigth = 1

# Create subplots
        fig, axes = plt.subplots(nrows=len(dd[0]), ncols=1, sharex=False, figsize=(8,6))
        plt.xlabel(x_label)
        # Plot each dataset
        for i, ax in enumerate(axes):
            ax.set_yticks([])
            if i != len(axes):
                ax.set_xticklabels([])
            if i == 0:
                ax.set_title(dd[1])
            d = dd[0][i]
            bp_resolution = bp_resolution_list[i]
            ax.barh(domain_bar_heigth,(d["end"]-d["start"]+1),left = (d["start"]-0.5),height = domain_bar_heigth,color=d["color"],align = "edge")
            if i ==  (len(axes)-1):
                my_step = round(tick_step/bp_resolution)
                tick_positions = np.arange(0,max(d["end"])+1,step=my_step)
                
                ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(tick_positions))
                ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(tick_labels))

        custom_legend = [Rectangle((0, 0), 1, 1, color='palevioletred', label='Overwound'),
                        Rectangle((0, 0), 1, 1, color='cornflowerblue', alpha=0.5, label='Underwound'),
                        Rectangle((0, 0), 1, 1, color='lightgrey', alpha=0.5, label='Neutral')]
# Create the legend
        plt.legend(custom_legend,["Overwound","Underwound","Neutral"],loc='upper left',bbox_to_anchor=(1.04, 1))
        plt.tight_layout()
        plt.savefig("Domains"+dd[1]+".png")
        #plt.show()
    
def get_x_tick_step(chrom_size):
    ##calculate numbers on x axis
    if chrom_size>1000 and chrom_size<=10000:
        tick_step = 1000
        x_label = "$10^3 bp$"
    if chrom_size>10000 and chrom_size<=100000:
        tick_step = 10000
        x_label = "$10^4 bp$"
    if chrom_size>100000 and chrom_size<=1000000:
        tick_step = 100000
        x_label = "$10^5 bp$"
    if chrom_size>1000000 and chrom_size<=10000000:
        tick_step = 1000000
        x_label = "$10^6 bp$"
    if chrom_size>10000000 and chrom_size<=100000000:
        tick_step = 10000000
        x_label = "$10^7 bp$"
    if chrom_size>100000000 and chrom_size<=1000000000:
        tick_step = 100000000
        x_label = "$10^8 bp$"
    if chrom_size<1000 or chrom_size>1000000000:
        print(chrom_size)
        raise Exception("Can't plot this size of chromosome!")
    return tick_step, x_label
 
def plot_signal(signal_list,domains_list,chrom_names,chrom_sizes,bp_resolution):
    signal_domains_list = zip(signal_list, domains_list,chrom_names,chrom_sizes)
    for s_d in signal_domains_list:
        sample_mean = np.mean(s_d[0])
        sample_std = np.std(s_d[0])
        
        sample_mean_plus = np.mean(np.take(s_d[0],np.where(s_d[0]>=np.float64(0))))
        sample_std_plus = np.std(np.take(s_d[0],np.where(s_d[0]>=np.float64(0))))
        sample_mean_minus = np.mean(np.take(s_d[0],np.where(s_d[0]<np.float64(0))))
        sample_std_minus = np.std(np.take(s_d[0],np.where(s_d[0]<np.float64(0))))

        fig, ax = plt.subplots()
        
    ##calculate numbers on x axis
        tick_step, x_label = get_x_tick_step(s_d[3])
        # define tick positions
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5.0))
        plt.xlabel(x_label)
        plt.ylabel("signal")
    
        tick_labels = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
        tick_positions = np.arange(0,len(s_d[0])+1,step=round(tick_step/bp_resolution))
        
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(tick_positions))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(tick_labels))

        ax.plot(s_d[0])
        ax.set_title(s_d[2])

        ax.bar(range(0,len(s_d[0])),s_d[0],width = 0.8)
        ax.axhline(y=sample_mean_plus+sample_std_plus,color="grey",linestyle = "--")
        ax.axhline(y=sample_mean_minus-sample_std_minus,color="grey",linestyle = "--")
        ax.axhline(y=0,color="black")
    ##Set domain bar height to 1/10 of the difference between max and min signal value
        domain_bar_heigth = (max(map(abs,s_d[0]))+min(map(abs,s_d[0])))/10
        ax.barh(min(s_d[0])-domain_bar_heigth,(s_d[1]["end"]-s_d[1]["start"]+1),left = (s_d[1]["start"]-0.5),height = domain_bar_heigth,color=s_d[1]["color"],align = "edge")
        plt.savefig("SignalAndDomains"+s_d[2]+".png")
        #plt.show()

def get_signal_bw(file_name_1,file_name_2, list_of_chrom,bp_resolution):
    bw_1 = pyBigWig.open(file_name_1)
    bw_2 = pyBigWig.open(file_name_2)
    list_of_signal = []
    chrom_names = []
    chrom_sizes = []

    nBasesCovered = bw_1.header()["nBasesCovered"]

    if not bw_1.isBigWig():
        print("Provided file 1 is not in BigWig format!")
    if not bw_2.isBigWig():
        print("Provided file 2 is not in BigWig format!")
    chroms_1 = bw_1.chroms()
    chroms_2 = bw_2.chroms()

    ##Make sure that provided files have the same number of chromosomes and chromosome length
    if chroms_1 != chroms_2:
        print("Provided files have different chromosome names and cannot be compared!")
        print(chroms_1)
        print(chroms_2)
        raise Exception("Bad input.")

    ##Check that bp resolution is at least half the chromosome length for all chromosomes
    if not list_of_chrom:
        for chrom_name, number_of_bp in chroms_1.items():
            ##Check that bp resolution is at least half the chromosome length for all chromosomes
            if number_of_bp < (bp_resolution*2):
                print(chrom_name)
                print(bp_resolution)
                print(number_of_bp)
                print("Requested resolution is too big for the chromosome size. Decrease bp resolution.")
                raise Exception("Bad input.")

            signal_1 = bw_1.stats(str(chrom_name),0,number_of_bp,nBins = round(number_of_bp/bp_resolution))
            signal_2 = bw_2.stats(str(chrom_name),0,number_of_bp,nBins = round(number_of_bp/bp_resolution))
            signal = [xi - yi for xi, yi in zip(signal_1, signal_2)]
            #print(signal)
            list_of_signal.append(signal)
            chrom_names.append(chrom_name)
            chrom_sizes.append(number_of_bp)
    else:
        for chrom_name in list_of_chrom:
            if chrom_name not in chroms_1.keys():
                print("This chromosome name is not among the chromosomes in the bw file.")
                print(chrom_name)
                raise Exception("Bad input.")
            else:
                number_of_bp = chroms_1[chrom_name]
                ##Check that bp resolution is at least half the chromosome length for all chromosomes
                if number_of_bp < (bp_resolution*2):
                    print(chrom_name)
                    print(bp_resolution)
                    print(number_of_bp)
                    print("Requested resolution is too big for the chromosome size. Decrease bp resolution.")
                    raise Exception("Bad input.")

                signal_1 = bw_1.stats(str(chrom_name),0,number_of_bp,nBins = round(number_of_bp/bp_resolution))
                signal_2 = bw_2.stats(str(chrom_name),0,number_of_bp,nBins = round(number_of_bp/bp_resolution))
                signal = [xi - yi for xi, yi in zip(signal_1, signal_2)]
                signal = list(map(lambda x:x*nBasesCovered, signal))
 
                list_of_signal.append(signal)
                chrom_names.append(chrom_name)
                chrom_sizes.append(number_of_bp)

    return list_of_signal, chrom_names, chrom_sizes

##Get signal and generate domains for each resolution
dd = []
for bp in bp_resolution_list:
    ll, chrom_names, chrom_sizes = get_signal_bw(file_1, file_2, list_of_chrom,bp)
    dd.append(find_domains(ll,bp))

##Create list with entry for each chromosome; each element in this list is a list of domains for each bp resolution
ddd = []
j=0
for c in chrom_names:
    cc=[] ## empty dummy list
    ddd.append(cc)
    for d in dd:
        ddd[j].append(d[j])
    j=j+1


if len(dd)==1:
    plot_signal(ll,dd[0],chrom_names,chrom_sizes,bp_resolution_list[0])
else:
    plot_domains_multiple_resolution(ddd,chrom_names,chrom_sizes,bp_resolution_list)

##Output domain coordinates into a table. Creates one file for each chromosome and resolution.
for i in range(0,len(chrom_names)):
    for j in range(0,len(bp_resolution_list)):
        ddd[i][j]['start'] = ddd[i][j]['start'] * bp_resolution_list[j]
        ddd[i][j]['end'] = ddd[i][j]['end'] * bp_resolution_list[j]
        ddd[i][j].drop(columns=['color'], inplace=True)
        ddd[i][j].to_csv(str(chrom_names[i])+"_"+str(bp_resolution_list[j])+"_domains.txt",sep='\t',index=False)
        

    