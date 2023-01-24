import roadrunner
rr = roadrunner.RoadRunner("https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000010.2?filename=BIOMD0000000010_url.xml")
results = rr.simulate(0, 2000, 200)
rr.plot()
print(results)