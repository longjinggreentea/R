###--- Import .csv or excel workbook data into R using {RSQLite} ---###


# https://www.r-bloggers.com/r-and-sqlite-part-1/
# 


library(RSQLite)
library(XLConnect)
library(sqldf)

if(file.exists('Test.sqlite')==TRUE) file.remove('Test.sqlite')
if(file.exists('Test1.sqlite')==TRUE) file.remove('Test1.sqlite')
if(file.exists('Tes2t.sqlite')==TRUE) file.remove('Test2.sqlite')

# create a data base
db<-dbConnect(SQLite(), dbname='Test.sqlite')

# sqldf("attach 'Test1.sqlite' as new")

# add data to the database - the hard way
dbSendQuery(conn = db, 
            "CREATE TABLE School
              ( SchID INTEGER,
                LOCATION TEXT,
                Authority TEXT,
                SchSIZE TEXT)")
dbSendQuery(conn=db,
            "INSERT INTO School
             VALUES (1, 'urban', 'state', 'medium')"
            )
dbSendQuery(conn=db,
            "INSERT INTO School
             VALUES (2, 'urban', 'independent', 'large')"
            )
dbSendQuery(conn=db,
            "INSERT INTO School
            VALUES (3, 'rural', 'state', 'small')"
            )
dbListTables(db)
dbListFields(db, "School")
dbReadTable(db,"School")
# It is a tedious way to enter data. Remove this table.
dbRemoveTable(db, "School")
dbListTables(db)

dbRemoveTable(db, "School")
dbRemoveTable(db, "Student")
dbRemoveTable(db, "Class")
dbListTables(db)

# Entering data from csv files
setwd("~/Data_study/R_and_SQLite")
list.files()

dbWriteTable(conn=db, name="Student", value="student.csv",
             row.names=FALSE, header=TRUE)
dbWriteTable(conn=db, name="Class", value="class.csv", 
             row.names=FALSE, header=TRUE)
dbWriteTable(conn=db, name="School", value="school.csv", 
             row.names=FALSE, header=TRUE)
dbReadTable(db, "Student"); dbListFields(db, "Student")
dbReadTable(db, "School")
dbReadTable(db, "Class")

dbRemoveTable(db, "School")

dbListTables(db)


# Import data frames into database

School <- read.csv("school.csv")

dbWriteTable(conn=db, name="School", value = School, row.name=FALSE)
dbReadTable(db, "School")

# Entering data from Excel workbooks
dbRemoveTable(db, "School")
dbRemoveTable(db, "Student")
dbRemoveTable(db, "Class")
dbListTables(db)

dbRemoveTable(db, "Student")
dbRemoveTable(db, "Class")
dbListTables(db)

# Import the three sheets of the Excel workbook
wb <- loadWorkbook("Test.xlsx")

Tables <- readWorksheet(wb, sheet = getSheets(wb))
str(Tables)
names(Tables)
dimnames(Tables$Student)

with(Tables, {
  dbWriteTable(conn = db, name = "Student", value = Student, row.names = FALSE)
  dbWriteTable(conn = db, name = "Class", value = Class, row.names = FALSE)
  dbWriteTable(conn = db, name = "School", value = School, row.names = FALSE)
})
dbListTables(db)
dbListFields(db,"Student")
dbReadTable(db,"Student")

dbDisconnect(db) # Close connection

rm(list = c("Tables", "Class", "School", "Student"))   # Remove data frames



