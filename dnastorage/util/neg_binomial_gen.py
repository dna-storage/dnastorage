#This file contains code for generating random numbers based on specified parameters, uses the inverse transform method with cutpoints to speed up generation
from scipy import special
from scipy import stats
import random
import math
import numpy




#Instantiate this class in order to generate random variables from the negative binomial distribution with gamma function parameterization (rather than classical parameterization)
class neg_bin:
    def __init__(self,mean,var):
        self._cumulative_array=[]
        self._prob_array=[]
        #variance should be greater than the mean for neg binomial, otherwise use Poisson distribution
        assert var>mean
        self._theta=mean**2/(var-mean)
        self._mean=mean
        self._var=var
        last_cumlative=0
        index=0
        #build the cummulative probability array and probability array 
        while last_cumlative < 99.999E-2:
            prob=self.pmf(index)
            self._prob_array.append(prob)
            if not self._cumulative_array:
                self._cumulative_array.append(prob)
            else:
                self._cumulative_array.append(self._cumulative_array[-1]+prob)
            if self._cumulative_array[-1]>=1:
                break
            index=index+1
            last_cumlative=self._cumulative_array[-1]
        if self._cumulative_array[-1]<1:
            self._cumulative_array[-1]=1

        #need to define cutoff points used for quickly getting random values
        self._m=len(self._cumulative_array)
        k=-1
        m_q=0
        self._cutoff_array=[]
        for j in range(0,self._m):
            while m_q<=j:
                k=k+1
                m_q=self._m*self._cumulative_array[k]
            self._cutoff_array.append(k)
            
        assert len(self._cumulative_array)==self._m
        
        
    #calculate a probability given a value
    def pmf(self,X):
        theta=self._theta
        mean=self._mean
        #gammaln is the natural log of the gamma function at the specified value, used to make sure the value from the gamma function does not explode
        log_prob=special.gammaln(X+theta)-special.gammaln(X+1)-special.gammaln(theta)+theta*(math.log(theta)-math.log(mean+theta))+X*(math.log(mean)-math.log(theta+mean))
        prob = math.exp(log_prob)
        return prob

    #function for generating negative binomial random value from a uniform random variable
    def gen(self):
        #Use this uniform random variable to 
        U=random.uniform(0,1)
        cumulative_index=self._cutoff_array[int(math.floor(U*self._m))]
        while U>self._cumulative_array[cumulative_index]:
            cumulative_index=cumulative_index+1
        return cumulative_index

    def get_mean(self):
        return self._mean

    def get_var(self):
        return self._var

    

def bins_array(data,bin_size):
    num_bins=int(math.ceil((max(data)-min(data))/bin_size))
    lower_bin=min(data)
    upper_bin=lower_bin+bin_size
    bin_array={}
    for bin_number in range(0,num_bins):
        frequency=0
        for count in data:
            if count>=lower_bin and count<upper_bin:
                frequency=frequency+1
        bin_array[lower_bin]=frequency
        lower_bin=upper_bin
        upper_bin=upper_bin+bin_size
    return bin_array,num_bins





def calculate_expected(num_samples,start,bin_size,prob_array,num_bins):
    expected={}
    index=start
    for i in range(0,num_bins):
        expected[int(index)]=num_samples*sum(prob_array[index:index+bin_size])
        index=index+bin_size
    return expected





    
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    _rand=neg_bin(200,3000)
    rand_reads=[]
    for i in range (0,1000000):
        rand_reads.append(int(_rand.gen()))

    #keep bin size integer, since we have a discrete random variable


    bin_size=1


    bins_,num_bins=bins_array(rand_reads,bin_size)
    start_point=min(rand_reads)
    probability_array=_rand._prob_array
    num_samples=len(rand_reads)
    expected=calculate_expected(num_samples,start_point,bin_size,probability_array,num_bins)
    assert len(expected)==len(bins_)
    sorted_bin_results=[bins_[x] for x in sorted(bins_)]
    sorted_expected_results=[expected[x] for x in sorted(expected)]
    chi_dist=0
    for index, value in enumerate(sorted_bin_results):
        #print "{}   {}".format(sorted_bin_results[index],sorted_expected_results[index])
        chi_dist=chi_dist+((value-sorted_expected_results[index])**2/(sorted_expected_results[index]))
        #print"chi {}".format(chi_dist)
    print (chi_dist)
    print (len(sorted_bin_results)-1-2)
    print (1 - stats.chi2.cdf(chi_dist, len(sorted_bin_results)-1-2))

    
    plt.hist(rand_reads,num_bins)
    plt.show()
