from DataBase import SQLdb


#execute = SQLdb

#execute.create_table('inulinases.db')

def clean_csv(file):
    list_strings = []
    with open(file, 'r') as file:
        lines = file.readlines()
        for i in lines:
            x = i.rstrip('\n').split(',')
            x = [item for item in x if item != '']
            x = ', '.join(x)
            list_strings.append(x)

    return list_strings


print(clean_csv('inulinases.csv'))

#execute.insert_dataset()
