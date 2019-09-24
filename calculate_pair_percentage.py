from itertools import permutations

protein_letters = 'ACDEFGHIKLMNPQRSTVWY'
combinations = list()
length = len(list(protein_letters))
for letter in range(length):
    for laterLetter in range(letter,length): #change to letter+1 to exclude pairs
        combinations.append(protein_letters[letter]+protein_letters[laterLetter])


print(combinations)
