import vcf 
import subprocess
import sys 
import os.path
import wave 
import tabix

#REQUIRED#
#pysam
#PyVcf
#pytabix


#myFile = sys.argv[1]
#myVCF = vcf.Reader(open(myFile, 'r'))

#myFile = sys.argv[1]
myFile = "clinvar.vcf"


print "Building BST(binary search tree) from "+myFile

print "This will take a few minutes "

# A utility class that represents an individual node in a BST 
class Node: 
	def __init__(self,key,info): 
		self.left = None
		self.right = None
		self.val = key
		self.info =  info 

# A utility function to insert a new node with the given key 
def insert(root,node): 
	if root is None: 
		root = node 
	else:
		while 1 == 1: 
			
			if root.val < node.val: 
				if root.right is None: 
					root.right = node
					return 
				else:
					root=root.right 
					#insert(root.right, node) 
			else:
				#print "left" 
				if root.left is None: 
					root.left = node 
					#print "new Node Left!"
					return
				else:
					root=root.left 
					#insert(root.left, node) 


# A utility function to search a given key in BST 
def search(root,key): 
     
     while 1 == 1: 
	    # Base Cases: root is null or key is present at root 
	    if root is None or root.val == key: 
	        return root 
	  
	    # Key is greater than root's key 
	    if root.val < key: 
	        root=root.right
	        #return search(root.right,key) 
	    
	    else:
	    	root=root.left
	    	# Key is smaller than root's key 
	    	#return search(root.left,key) 

    # Python program to demonstrate insert operation in binary search tree 

# A utility function to do inorder tree traversal 
def inorderLeft(root): 

	while root != None: 
		#inorder(root.left) 
		print(root.val) 
		#inorder(root.right) 
		root=root.left

def inorderRight(root): 
	

	while root != None: 
		#inorder(root.left) 
		print(root.val) 
		#inorder(root.right) 
		root=root.right

bigDict={}
def buildDataStore():
	#sys.setrecursionlimit(8000)
	lastChrom=0
	with open(myFile) as g:
			lineCount=0
			
			for line in g:

					lineCount +=1
					if lineCount<29:
						continue
					

					if lineCount%10000==0:
						print "lines processed: "+str(lineCount)

					myList=line.split("\t")
					if len(myList)>1:
						if myList[0] > lastChrom:
							print "Processing Chromosome: "+myList[0]
							lastChrom=myList[0]
						r=Node(125000000,[])

						if myList[0] in bigDict:
							r=bigDict[myList[0]]
						
						myLoc=int(myList[1])
						insert(r,Node(myLoc,myList))
						bigDict[myList[0]]=r



#main test function that returns boolean for whether disease was found at location, and a list of the associated features 

def testForDisease(chromosome, position):
	myResult=search(bigDict[str(chromosome)], position)
	if myResult == None:
		return False , []
	return True, myResult.info


buildDataStore()

# Driver program to test the above functions 
# Let us create the following BST 
#	 50 
# /	 \ 
# 30	 70 
# / \ / \ 
# 20 40 60 80 
'''
r = Node(50) 
insert(r,Node(30)) 
insert(r,Node(20)) 
insert(r,Node(40)) 
insert(r,Node(70)) 
insert(r,Node(60)) 
insert(r,Node(80)) 

# Print inoder traversal of the BST 
inorder(r) 

# This code is contributed by Bhavya Jain 
'''
