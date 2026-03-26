from scipy.special import roots_legendre

n_total = 100

with open("gl_to_cpp.txt", "w") as file:
    for n in range(2, n_total + 1):
        roots, weights = roots_legendre(n)
        print(f"case {n}: " + "{", file=file)
        print("weights_ = {", file=file)
        for weight in weights[:-1]:
            print(f"{weight}, ", file=file)
        print(f"{weights[-1]}" + "};", file=file)

        print("abscissae_ = {", file=file)
        for root in roots[:-1]:
            print(f"{root}, ", file=file)
        print(f"{roots[-1]}" + "};", file=file)
        print("};", file=file)
