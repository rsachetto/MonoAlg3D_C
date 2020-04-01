import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_activation_time_results (data_sc1,data_sc2):
    
    labels = np.arange(1,21)    # Label name
    x = np.arange(1,21)         # Label position
    width = 0.35                # Bar width

    #fig, ax = plt.subplots()
    #ax.bar(x - width/2,data_sc1,width,label='SC1')
    #ax.bar(x + width/2,data_sc2,width,label='SC2')
    plt.bar(x - width/2,data_sc1,width,label='SC1')
    plt.bar(x + width/2,data_sc2,width,label='SC2')

    #ax.set_ylabel("RMS",fontsize=15)
    #ax.set_title("Error activation time",fontsize=15)
    #ax.set_xticks(x)
    #ax.set_xticklabels(labels)
    #ax.legend(loc=0,fontsize=10)
    plt.ylabel("RMS",fontsize=15)
    plt.title("Error activation time",fontsize=15)
    plt.xticks(x,labels)
    #plt.xticklabels(labels)
    plt.legend(loc=0,fontsize=10)

    #plt.show()
    #plt.savefig("error_activation_time.pdf")
    plt.savefig("error_activation_time.png",dpi=400)

def plot_apd_results (data_sc1,data_sc2):
    
    labels = np.arange(1,21)    # Label name
    x = np.arange(1,21)         # Label position
    width = 0.35                # Bar width

    #fig, ax = plt.subplots()
    #ax = fig.add_axes([0,0,1,1])
    #ax.bar(x - width/2,data_sc1,width,label='SC1')
    #ax.bar(x + width/2,data_sc2,width,label='SC2')
    plt.bar(x - width/2,data_sc1,width,label='SC1')
    plt.bar(x + width/2,data_sc2,width,label='SC2')

    #ax.set_ylabel("RMS",fontsize=15)
    #ax.set_title("Error APD",fontsize=15)
    #ax.set_xticks(x)
    #ax.set_xticklabels(labels)
    #ax.legend(loc=0,fontsize=10)
    plt.ylabel("RMS",fontsize=10)
    plt.title("Error APD",fontsize=15)
    plt.xticks(x,labels)
    #plt.xticklabels(labels)
    plt.legend(loc=0,fontsize=10)

    #plt.show()
    #plt.savefig("error_apd.pdf")
    plt.savefig("error_apd.png",dpi=400)

def main():
	
    if len(sys.argv) != 5:
        print("------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <result_filename>" % sys.argv[0])
        print("------------------------------------------------------------------------------------------------------")
        print("Example:> python %s ../outputs/elnaz-purkinje-coupled-errors/results-error-activation-map-sc1.txt" % (sys.argv[0]))
        print("                    ../outputs/elnaz-purkinje-coupled-errors/results-error-activation-map-sc2.txt")
        print("                    ../outputs/elnaz-purkinje-coupled-errors/results-error-apd-map-sc1.txt")
        print("                    ../outputs/elnaz-purkinje-coupled-errors/results-error-apd-map-sc2.txt")
        print("------------------------------------------------------------------------------------------------------")
        return 1

    error_at_sc1_filename = sys.argv[1]
    error_at_sc2_filename = sys.argv[2]
    error_apd_sc1_filename = sys.argv[3]
    error_apd_sc2_filename = sys.argv[4]

    error_at_sc1 = np.genfromtxt(error_at_sc1_filename)
    error_at_sc2 = np.genfromtxt(error_at_sc2_filename)
    error_apd_sc1 = np.genfromtxt(error_apd_sc1_filename)
    error_apd_sc2 = np.genfromtxt(error_apd_sc2_filename)

    plot_activation_time_results(error_at_sc1,error_at_sc2)
    plt.clf()
    plot_apd_results(error_apd_sc1,error_apd_sc2)

if __name__ == "__main__":
	main()
