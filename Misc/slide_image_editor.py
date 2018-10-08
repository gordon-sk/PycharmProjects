from PIL import Image
filename = 'picture1.png'
img = Image.open(filename)
pix = img.load()
print(pix[0,0])
print(img.size)

for x in range(img.size[0]):
    for y in range (img.size[1]):
        if pix[x,y] == (255, 255, 255, 255):
            pix[x,y] = (255, 255, 255, 0)

img.save("uhhh.png")