import albumentations as A
import os
import cv2
import matplotlib.pyplot as plt

SavePathImgs = r'augmented_data\imgs'
SavePathMasks = r'augmented_data\masks'
os.makedirs(SavePathImgs, exist_ok=True)
os.makedirs(SavePathMasks, exist_ok=True)

PathUsualDataImgs = r'original_annotated_data\imgs'
PathUsualDataLabels = r'original_annotated_data\labels'

ImgsList = os.listdir(PathUsualDataImgs)
LabelsList = os.listdir(PathUsualDataLabels)


def Transforms(p, sets):
    t = [A.RandomRotate90(p=0.8),
         A.Flip(p=0.8),
         A.Transpose(p=0.8),
         A.OneOf([
             A.IAAAdditiveGaussianNoise(),
             A.GaussNoise(),
         ], p=0.75),
         A.OneOf([
             A.MotionBlur(p=0.33),
             A.MedianBlur(blur_limit=3, p=0.33),
             A.Blur(blur_limit=3, p=0.33),
         ], p=0.75),
         A.ShiftScaleRotate(shift_limit=0.1, scale_limit=0.3, rotate_limit=45, p=0.5),
         A.OneOf([
             A.OpticalDistortion(p=0.3),
             A.GridDistortion(p=0.3),
             A.IAAPiecewiseAffine(p=0.3),
         ], p=0.75),
         A.OneOf([
             A.CLAHE(clip_limit=2),
             A.IAASharpen(),
             A.IAAEmboss(),
             A.RandomBrightnessContrast(),
         ], p=0.75),
         A.HueSaturationValue(p=0.5)]
    if not sets:
        t.append(A.OneOf([
            A.RandomCrop(256, 256),
        ], p=1))
    else:
        t.append(A.Resize(512, 512, always_apply=True))
        t.append(A.OneOf([
            A.RandomCrop(256, 256),
        ], p=1))
    return A.Compose(t, p=p)


for img, label in zip(ImgsList, LabelsList):
    Image = cv2.imread(os.path.join(PathUsualDataImgs, img))
    Label = cv2.imread(os.path.join(PathUsualDataLabels, label))
    data = {'image': Image, 'mask': Label}
    for i in range(800):
        Albu = Transforms(p=1, sets=False)
        Augmented = Albu(**data)
        cv2.imwrite(os.path.join(SavePathImgs, f'{img}_{i}.png'), Augmented['image'])
        cv2.imwrite(os.path.join(SavePathMasks, f'{label}_{i}.png'), Augmented['mask'])
