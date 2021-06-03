import sqlite3
from sqlite3 import Error

database = r"chembl_28/chembl_28_sqlite/chembl_28.db"

# create a database connection
try:
    conn = sqlite3.connect(database)
    with conn:
        cur = conn.cursor()
        cur.execute("SELECT * FROM drug_mechanism")
        rows = cur.fetchall()

        i=0
        for row in rows:
            i+=1
            print(row)
            if (i>10):
                break

except Error as e:
    print(e)