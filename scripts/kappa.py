import numpy as np

def calc_kappa (length,delta):
	return (pow(length,4)-pow(delta,4)) / (12.0*pow(length,2))

def main ():
	cell_length = 20.0
	delta = 30.0
	print("Kappa = %.10lf" % (calc_kappa(cell_length,delta)))

if __name__ == "__main__":
	main()
