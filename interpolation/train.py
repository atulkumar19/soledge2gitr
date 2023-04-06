import torch
import nn_learner

DB_PATH='SOLEDGE.db'

model =nn_learner.retrain(db_path=DB_PATH)

#Ensure model can be pickled and unpickled.
torch.save(model,"data.pt")

