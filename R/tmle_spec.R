#tmle_spec
#parameter(s)
#submodel + loss function(learner)
#dag (default)

#should be able to fit by specifying a minimum of things
#tmle_task
#sl3 learners for relevant factors (default to glm)

#node info is currently somewhat split across tmle_task and likelihood
#might be better off in just tmle_task
#then, can assume a default before a task is even specified
#how will this affect missingness?
#tmle_task needs this info to make regression tasks

#todo: make sure nonoutcome nodes get formatted correctly !!!important!!!
# (e.g. multinomial A should not be treated as continuous)
