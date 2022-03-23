# Step 0. 
# a. read the documentation of optim() function to 
#    understand it

# Step 1.
# a. create function to compute gradient

# Step 2.
# a. use BFGS quasi-Newton option of optim() function
#    to compute MLE of mu and sigma by
#    passing computed gradient into optim() function

# Step 3.
# a. generate 200 data points using the gen() provided
#    with mu = (-1,1,2)^T, sigma = (1,.7,.7,.7,1,.7,.7,.7,1)
# b. set seed to 2022


# Step 4.
# a. use initial values: sigma = ID matrix, mu = (0,0,0)^T
# b. abstol = 10^-5
# c. set an appropriate value for trace of optim()


#Let's all use tolerr= 10e-6 and tolgrad=10e-5.









