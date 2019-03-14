
import vcf 
import subprocess
import sys 
import os.path
import wave 
import tabix
from random import *
import getVariants2 as getVariants2
import pandas as pd


import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
import numpy as np


#REQUIRED#
#pysam
chromMax={1:248979990,2:242193529,3:198295559,4:190214555,5:181538259,
6:170805979,7:159345973,8:145138636,9:138394717,10:133797422,11:135086622,
12:133275309,13:114364328,14:107043718,15:101991189,16:90338345,17:83257441,
18:80373285,19:58617616,20:64444167,21:46709983,22:50818468}


#for i in [1,2,3,4,5,6,7,8]:
print chromMax[randint(1, 22)] 

def buildPseudoRandomTrainingFrame(chrom,start,window,size,spacing):
	myCols=["chrom","pos","window","disease","details"]
	theDF = pd.DataFrame(columns=myCols) #creates a new dataframe that's empty
	myChrom=chrom
	myArray=np.zeros(123,)
	print myArray.shape

	currentPostion=start
	for i in range(size):
		
		currentPostion=currentPostion+spacing
		#myLabel, myDisease=mutationData.testForDisease(chrom,newStart)

		#myArray=np.append(myArray,myLabel,axis=0)
		newArray=np.array([chrom,currentPostion,window])
		myVector=getVariants2.getFeatureVector(chrom, currentPostion, window)
		newArray=np.append(newArray,myVector)
		print newArray.shape
		print newArray
		myArray=np.append([myArray],[newArray],axis=0)

	return myArray

	#theDF = theDF.append(theDF, ignore_index = True)
myTrainFrame=buildPseudoRandomTrainingFrame(1,20000,10,100,5)
print(myTrainFrame)

def trainModel():
	print "yay"



window=100

#this will execute the mutationData.py file, which will build a BST from the clinvar.vcf file. Both files must be in the same directory. 
#this will take 3-10 minutes to build. 

myVector = getVariants2.getFeatureVector(4, 505000, window)
#import mutationData as mutationData

print myVector
#getVariants(4 ,505000, 10000)
#example: code below tests position 86952254 on chromosome 11 for a known diesease caused by mutation. Returns two objects, a boolean and a list of features related to the disease. 
#diseaseFound , diseaseInfo = mutationData.testForDisease(4,505000)
#print diseaseFound
#print diseaseInfo



batch_size = 2
num_classes = 2
epochs = 12

# input image dimensions
img_rows, img_cols = 10, 24

# the data, split between train and test sets
#(x_train, y_train), (x_test, y_test) = mnist.load_data()

#x_train = x_train.reshape(60000,img_cols,img_rows,1)
#x_test = x_test.reshape(10000,img_cols,img_rows,1)

print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print(x_test.shape[0], 'test samples')

# convert class vectors to binary class matrices
y_train = keras.utils.to_categorical(y_train, num_classes)
y_test = keras.utils.to_categorical(y_test, num_classes)
'''
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
                 activation='relu',
                 input_shape=(28,28,1)))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))

model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])

model.fit(x_train, y_train,
          batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test))
score = model.evaluate(x_test, y_test, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
'''
