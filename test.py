

def recur(i,fin,out = 0):
    if i > fin :
        return (out)
    out += i + 1
    return recur(i+1 , fin , out)

print(recur(1,10))