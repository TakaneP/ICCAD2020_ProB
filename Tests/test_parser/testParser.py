#!/usr/bin/python3
import subprocess

def main():
    with open("test", "rt") as f:
        for line in f:
            line = line.rstrip()
            test = line.split(" ")
            x = test[0]
            y = test[1]
            z = test[2]
            demand = int(test[3])
            
            out = subprocess.check_output(["./evaluator/evaluator","cases/case1.txt","--dbgrow",x,"--dbgcol",y,"--dbglay",z])
            out = out.decode("utf8")
            golden = 0
            if "current total dmd:" in out:
                out = out.split(":")[-2]
                golden = int(out.split("T")[0])
            if golden != demand:
                print("Demand Wrong at " + x + " " + y + " " + z)
    




if __name__ == "__main__":
    main()

