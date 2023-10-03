from DataBase import SQLdb


#execute = SQLdb

#execute.create_table('inulinases.db')

with open('inulinases.csv', 'r') as file:
    lines = file.readlines()
    for i in lines:
        x = i.rstrip('\n').split(',')
        if '' in x:
            print(x)
#execute.insert_dataset()
