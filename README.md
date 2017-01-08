# ParticleJetClassification
Individual Project where I combined data science techniques to tackle a well known physics problem

Particle jets are the resulting spray of particles when quarks and gluons are made in the Large Hadron Collider. In order to learn about the fundamental particles, we are forced to learn from the jet that they produce. Unfortunately, discriminating bewtween different types of jets is a non-trivial process. I seek to add to the volume of work on this problem by applying various machine learning techniques to attempt to classify the jets. 

The simulated data used for this project was generated using the Madgraph5 software, which produced the information of the initial fundamental particles after a collision. Then, I used pythia to shower these initial particles, and FastJet to group the shower into identifiable jets. I then matched the initial particles to the jets based on angular position to get the "correct" classification. 

Then, I was able to perform some analysis and classification in R using supervised learning techniques. The results were promising for 2/3 jet classifications.
