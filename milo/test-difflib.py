from difflib import *;
    
#      0    5    10   15   20
ref = "AAATTAAAcatCATCATTAGTAGTACCGGGG"
chk = "AAGCTAAAGcatTTTTTTAGTAGCGGGG"
s = SequenceMatcher(None, ref, chk, False)
matches = s.get_matching_blocks()

print("Ref: " + ref)
print("Chk: " + chk)
finSequence = []
lastM = None
aIndex = 0
bIndex = 0
for m in matches:
    matchedString = ref[m.a:m.a+m.size]
    finSequence.append({'type': 'm', 'str': matchedString})

    if lastM == None: # First match
        aIndex = m.a + m.size
        bIndex = m.b + m.size
        lastM = m;
        continue;
    
    aSize = m.a - aIndex
    bSize = m.b - bIndex
    
    if aSize == bSize:
        seq = chk[bIndex:bIndex+bSize]
        print("Sub: " + seq)
    elif aSize == 0 and bSize > 0:
        seq = chk[bIndex:bIndex+bSize]
        print("Ins: " + seq)
    elif bSize == 0 and aSize > 0:
        seq = ref[aIndex:aIndex+aSize]
        print("Del: " + seq)
    elif bSize > aSize:
        seq = chk[bIndex:bIndex+bSize]
        print("InsSub: " + seq)
    elif aSize > bSize:
        seq = chk[bIndex:bIndex+bSize]
        print("DelSub: " + seq)
        
    print ('aP: ' + str(aIndex) + ' bP: ' + str(bIndex))
    print ('aD: ' + str(aSize) + ' bD: ' + str(bSize))
    
    
    print (matchedString)
    print (m)
    aIndex = m.a + m.size
    bIndex = m.b + m.size
    lastM = m;

print ( finSequence )