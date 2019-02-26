from subsbml import *

model = createBasicSubsystem('cell')
m = SimpleModel(model.getSBMLDocument().getModel())
m.createNewUnitDefinition('per_second',libsbml.UNIT_KIND_SECOND,-1)
m.createNewUnitDefinition('count',libsbml.UNIT_KIND_DIMENSIONLESS,-1)
m.createNewSpecies(['x','y','c1'],'cell',[10,100,0], False, 'count')
m.createNewParameter(['kc','k1'],[6,1], True, 'per_second')
m.createSimpleReaction('r1','x + y --> c1','kc*x*y')
m.createSimpleReaction('r2','c1 --> ','k1*c1')

libsbml.writeSBML(model.getSBMLDocument(), 'sbml_example.xml')