#########################
#Load R packages
#########################
library(lubridate)


############################
#Basic R Objects
############################

#Use "<-" to create any object 
x <- 5	 #create a scalar object x with value 5
x

#Use c() to concatenate entries into a vector
y <- c(1,2,3,4)  #create a numeric vector y = (1 2 3 4)
y
st <- c("NY", "CA", "FL", "TX") #vector st with strings
st

#Use matrix() to turn vector y into a 2x2 matrix z
z <- matrix(y, ncol = 2, nrow = 2) 
z

#list combines "y" and "st" (different data types)
#label columns "order" and "states" / change with names()
l <- list(order = y, states = st)
l

#dataframe combines "y" and "st" (different data types)
#label columns "order" and "states" / change with names()
df <- data.frame(order = y, states = st)
df

#access elements of vectors, matrices and data frames
y[2] #element 2 of vector "y"
z[1,2] #element in row 1 and column 2 of matrix "z"
df$states #vector "state" in "df" (same for lists)
df$states[1:3] #elements 1 to 3 in vector "state" of "df"
l$order[c(1,4)] #elements 1 and 4 in vector "order" of "l"

#check data type of objects
word <- as.character("test")
word
is.character(word)

number <- as.numeric(5)
number
is.numeric(number)

categorical <- as.factor(c("male","female"))
categorical
is.factor(categorical)

day <- as.Date("13/04/81", "%d/%m/%y")
is.Date(day)

response <- TRUE
is.logical(response)




############################
#Basic R operations
############################

head(df, n = 2) #show first 2 rows of data frame df
tail(df, n = 2) #show last 2 rows of data frame df

x - 2 #subtract scalar 2 from scalar x
y + y #add vector y to vector y (need same dimensions)
z - 3 #subtract scalar 3 from matrix z

x * 2 #multiply scalar 2 with scalar x
y * y #multiply vector y with itself (element by element!)
z * 3 #multiply scalar 3 with matrix z
z * z #multiply matrix z with itself (element by element!)

z %*% z #multiply matrix z with itself (matrix product)
y %*% t(y) #multiply vector y with itself (vector product)
#Required: no. columns of 1st = no. rows of 2nd matrix

det(z) #compute determinant of matrix z
diag(z) #return main diagonal of matrix z
diag(y %*% t(y))

solve(z) #invert matrix z
#Required: nonsingular matrix (determinant ??? 0)
t(z) #transpose matrix z


