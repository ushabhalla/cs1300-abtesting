from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2

from abtest_data import *

# You can comment out these lines! They are just here to help follow along to the tutorial.
# print(t_dist.cdf(-2, 20)) # should print .02963
# print(t_dist.cdf(2, 20)) # positive t-score (bad), should print .97036 (= 1 - .2963)

# print(chi2.cdf(23.6, 12)) # prints 0.976
# print(1 - chi2.cdf(23.6, 12)) # prints 1 - 0.976 = 0.023 (yay!)

# TODO: Fill in the following functions! Be sure to delete "pass" when you want to use/run a function!
# NOTE: You should not be using any outside libraries or functions other than the simple operators (+, **, etc)
# and the specifically mentioned functions (i.e. round, cdf functions...)

def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append

def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    #TODO: fill me in!
    s = 0
    for num in nums:
        s += num

    return s/len(nums)
# print("get_avg:")
# print(get_avg([0,1,2,3,4,5,6]))

def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    #TODO: fill me in!
    s = 0

    avg = get_avg(nums)
    for num in nums:
        s += (num-avg)**2

    return (s/(len(nums)-1))**(1/2)

def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    #TODO: fill me in!
    l_a = len(a)
    f = (get_stdev(a)**2)/l_a

    l_b = len(b)
    s = (get_stdev(b)**2)/l_b
    return (f + s)**(1/2)

def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    #TODO: fill me in!
    se = get_standard_error(a, b)**4
    l_a = len(a)
    f = (((get_stdev(a)**2)/l_a)**2)/(l_a-1)

    l_b = len(b)
    s = (((get_stdev(b)**2)/l_b)**2)/(l_b-1)

    return round(se/(f+s))

def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    #TODO: fill me in!
    ts =(get_avg(a)-get_avg(b))/get_standard_error(a,b)
    # if ts > 0:
    #     ts *= -1 
    return ts

def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    #TODO: fill me in!
    ts = get_t_score(a,b)
    if ts > 0:
        ts *= -1
    return t_dist.cdf(ts, get_2_sample_df(a,b))


# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().
# def row_sum(observed_grid, ele_row):
# def col_sum(observed_grid, ele_col):
# def total_sum(observed_grid):
# def calculate_expected(row_sum, col_sum, tot_sum):

def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    #TODO: fill me in!
    nrows = len(observed_grid)
    ncols = len(observed_grid[0])
    row_sums = [0]*nrows
    col_sums = [0]*ncols
    total = 0
    j = 0
    for r in range(0,nrows):
        i = 0
        for val in observed_grid[r]:
            row_sums[j] += val
            col_sums[i] += val
            i += 1
            total += val
        j += 1

    expected_grid = []

    for i in range(0,nrows):
        row = []
        for j in range(0,ncols):
            row.append((row_sums[i]*col_sums[j])/total)
        expected_grid.append(row)
    return expected_grid

def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    #TODO: fill me in!
    df = (len(observed_grid)-1)*(len(observed_grid[0])-1)
    if df < 0:
        df *= -1
    return df

def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    chi2 = 0
    expected_grid = get_expected_grid(observed_grid)
    for i in range(len(observed_grid)):
        for j in range(len(observed_grid[0])):
            chi2 += ((observed_grid[i][j] - expected_grid[i][j])**2)/expected_grid[i][j]
    return chi2

def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    #TODO: fill me in!
    return 1 - chi2.cdf(chi2_value(observed_grid), df_chi2(observed_grid))

# These commented out lines are for testing your main functions. 
# Please uncomment them when finished with your implementation and confirm you get the same values :)
def data_to_num_list(s):
  '''
    Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums. 
    This will be useful when you need to run your tests on your cleaned log data!
    :param str: string holding data
    :return: the spliced list of numbers
    '''
  return list(map(float, s.split()))


# Time to completion t test
a_ttc_list = data_to_num_list(a_time_to_completion)
b_ttc_list = data_to_num_list(b_time_to_completion)
print("Time to completion t test:")
print(get_t_score(a_ttc_list,b_ttc_list))
print(perform_2_sample_t_test(a_ttc_list, b_ttc_list))

# Return rate chi2 test
a_c1_list = data_to_num_list(a_count) 
b_c1_list = data_to_num_list(b_count)
c1_observed_grid = [a_c1_list, b_c1_list]
print("Return reate chi2 test:")
print(chi2_value(c1_observed_grid))
print(perform_chi2_homogeneity_test(c1_observed_grid))

# # t_test 1:
# a_t1_list = data_to_num_list(a1) 
# b_t1_list = data_to_num_list(b1)
# print("should be -129.500:")
# print(get_t_score(a_t1_list, b_t1_list)) # this should be -129.500
# print("should be 0.00:")
# print(perform_2_sample_t_test(a_t1_list, b_t1_list)) # this should be 0.0000
# # why do you think this is? Take a peek at a1 and b1 in abtesting_test.py :)

# # t_test 2:
# a_t2_list = data_to_num_list(a2) 
# b_t2_list = data_to_num_list(b2)
# print("should be -1.48834:")
# print(get_t_score(a_t2_list, b_t2_list)) # this should be -1.48834
# print("should be .082379:")
# print(perform_2_sample_t_test(a_t2_list, b_t2_list)) # this should be .082379

# # t_test 3:
# a_t3_list = data_to_num_list(a3) 
# b_t3_list = data_to_num_list(b3)
# print("should be -2.88969:")
# print(get_t_score(a_t3_list, b_t3_list)) # this should be -2.88969
# print("should be .005091:")
# print(perform_2_sample_t_test(a_t3_list, b_t3_list)) # this should be .005091
# """
# """
# # chi2_test 1:
# a_c1_list = data_to_num_list(a_count_1) 
# b_c1_list = data_to_num_list(b_count_1)
# c1_observed_grid = [a_c1_list, b_c1_list]
# print("should be 4.103536:")
# print(chi2_value(c1_observed_grid)) # this should be 4.103536
# print("should be .0427939:")
# print(perform_chi2_homogeneity_test(c1_observed_grid)) # this should be .0427939

# # chi2_test 2:
# a_c2_list = data_to_num_list(a_count_2) 
# b_c2_list = data_to_num_list(b_count_2)
# c2_observed_grid = [a_c2_list, b_c2_list]
# print("should be 33.86444:")
# print(chi2_value(c2_observed_grid)) # this should be 33.86444
# print("should be 0.0000:")
# print(perform_chi2_homogeneity_test(c2_observed_grid)) # this should be 0.0000
# # Again, why do you think this is? Take a peek at a_count_2 and b_count_2 in abtesting_test.py :)

# # chi2_test 3:
# a_c3_list = data_to_num_list(a_count_3) 
# b_c3_list = data_to_num_list(b_count_3)
# c3_observed_grid = [a_c3_list, b_c3_list]
# print("should be .3119402:")
# print(chi2_value(c3_observed_grid)) # this should be .3119402
# print("should be .57649202:")
# print(perform_chi2_homogeneity_test(c3_observed_grid)) # this should be .57649202