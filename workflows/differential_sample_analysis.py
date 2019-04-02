from software.differentialexpression import DifferentialExpression
samples = ["Y7640", "Y7652", "Y7668", "Y8841"]


print ("Running DE")
de = DifferentialExpression(samples[:2])
de.run()
