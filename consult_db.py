from DataBase import SQLdb


if __name__ == '__main__':
    execute = SQLdb()
    execute.consult_select()
    column_index = 3  # Índice da terceira coluna (começando em 0)
    values = execute.select_all_by_column_index(column_index)
    for i in values:
        print(i)
