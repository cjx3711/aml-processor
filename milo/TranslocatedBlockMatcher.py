from difflib import SequenceMatcher
from DataTypes import *

class TranslocatedBlockMatcher:
	def __init__(self):
		pass
		
	def findLongestMatch(self, partitionsOnRef, partitionsOnRead, availableOnRead, seqMatcher):
		matchList = []
		for i in range(len(partitionsOnRead)):
			minPartOnRead, maxPartOnRead = partitionsOnRead[i]
			matchListPerPart = [] # List of (low, high) coordinates available for each partition to be matched against
			for minAvailable, maxAvailable in availableOnRead:
				if minAvailable > maxPartOnRead: # If current segment is entirely beyond partition bounds, break
					break
				elif maxAvailable < minPartOnRead: # If current segment is entirely before partition bounds, continue
					continue
				else: # If current segment is entirely within partition bounds, append and continue
					matchListPerPart.append((minAvailable, maxAvailable))

			minPartOnRef, maxPartOnRef = partitionsOnRef[i]
			if matchListPerPart: # Find longest common substring between a partition and all applicable locations on read
				matchList.append(max([seqMatcher.find_longest_match(readLow, readHigh, minPartOnRef, maxPartOnRef) for readLow, readHigh in matchListPerPart], key = lambda x: x.size))
			else: # If there are no applicable locations on read, append an EmptyMatch
				matchList.append(EmptyMatch(0, 0, 0))
		# Return the longest common substring across all partitions
		return max(matchList, key = lambda x: x.size)

	def exciseAvailable(self, excisionStart, excisionSize, allRegions):
		regionsAftExcision = []
		excisionEnd = excisionStart + excisionSize
		for reg in allRegions:
			if excisionStart > reg[1] or excisionEnd < reg[0]: # If we are not excising from current region, append region unaltered
				regionsAftExcision.append(reg)
			else:
				if reg[0] < excisionStart: # If we fail to excise everything from the front of current region, create new region between the front and excisionStart
					regionsAftExcision.append((reg[0], excisionStart))
				if reg[1] > excisionEnd: # If we fail to excise everything from the back...
					regionsAftExcision.append((excisionEnd, reg[1]))
				# If we excise the entire region, nothing is appended
		return regionsAftExcision

	def excisePartitions(self, excisionData, partitionsOnRead, partitionsOnRef):
		regsAftExcisionRead = []
		regsAftExcisionRef = []
		excisionEndRead = excisionData.a + excisionData.size
		excisionEndRef = excisionData.b + excisionData.size

		for regRead, regRef in zip(partitionsOnRead, partitionsOnRef):
			if excisionData.a > regRead[1] or excisionEndRead < regRead[0]: # If we are not excising from current region, append region unaltered
				regsAftExcisionRead.append(regRead)
				regsAftExcisionRef.append(regRef)
			else:
				# When excising, iff we do not completely excise every base from the ends of both reference and read partitions, a new, smaller partition is formed at said end
				if regRead[0] < excisionData.a and regRef[0] < excisionData.b:
					regsAftExcisionRead.append((regRead[0], excisionData.a))
					regsAftExcisionRef.append((regRef[0], excisionData.b))
				if regRead[1] > excisionEndRead and regRef[1] > excisionEndRef:
					regsAftExcisionRead.append((excisionEndRead, regRead[1]))
					regsAftExcisionRef.append((excisionEndRef, regRef[1]))
				# If every base is completely excised from the end of either reference or read partition, or both, no new partition is added
		return regsAftExcisionRead, regsAftExcisionRef

	def findTranslocatedMatchingBlocks(self, read, ref1, ref2):
		# Sets the coordinates of the partitions of both references onto the read
		partitionsOnRefR1 = [(0, len(ref1))]
		partitionsOnReadR1 = [(0, len(read))]
		seqMatcherR1 = SequenceMatcher(None, read, ref1, False)

		partitionsOnRefR2 = [(0, len(ref2))]
		partitionsOnReadR2 = [(0, len(read))]
		seqMatcherR2 = SequenceMatcher(None, read, ref2, False)

		availableOnRead = partitionsOnReadR1[:] # Bases available for matching shared between both reads (i.e. the same part of a read cannot be assigned to parts of both references)
		matchingBlocks = []

		while True:
			# Calculate the longest match between the read and both references
			longestMatchR1 = self.findLongestMatch(partitionsOnRefR1, partitionsOnReadR1, availableOnRead, seqMatcherR1)
			longestMatchR2 = self.findLongestMatch(partitionsOnRefR2, partitionsOnReadR2, availableOnRead, seqMatcherR2)

			if longestMatchR1.size == 0 and longestMatchR2.size == 0: # If we no longer have a single matching block, break and return all previous matching blocks
				break
			elif longestMatchR1.size > longestMatchR2.size: # Pick the better matching block and confirm it as the match for this iteration, excise the block from the partitions and available bases
				partitionsOnReadR1, partitionsOnRefR1 = self.excisePartitions(longestMatchR1, partitionsOnReadR1, partitionsOnRefR1)
				matchingBlocks.append(("R1", longestMatchR1))
				availableOnRead = self.exciseAvailable(longestMatchR1.a, longestMatchR1.size, availableOnRead)
			else:
				partitionsOnReadR2, partitionsOnRefR2 = self.excisePartitions(longestMatchR2, partitionsOnReadR2, partitionsOnRefR2)
				matchingBlocks.append(("R2", longestMatchR2))
				availableOnRead = self.exciseAvailable(longestMatchR2.a, longestMatchR2.size, availableOnRead)

			for emptyChecker in [partitionsOnReadR1, partitionsOnRefR1, partitionsOnReadR2, partitionsOnRefR2, availableOnRead]: # Put EmptyMatch objects into empty partitions, so even when one reference has exhausted, the other can continue matching
				if not emptyChecker:
					emptyChecker.append((0,0))
		return matchingBlocks