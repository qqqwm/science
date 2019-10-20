import math
import matplotlib.pyplot as plt
m_N = 0.94
def energy_finder(p_0, n = 2, m_N = 0.94):
    b = -2*m_N*(n*p_0**2 + (n**2+1)*m_N*math.sqrt(m_N**2 + p_0**2) + 2*n*m_N**2)
    a = (n**2+1)*m_N**2 + 2*n*m_N*math.sqrt(m_N**2 + p_0**2)
    c = (n**2+1)*m_N**4 + (n**2+1)*m_N**2*p_0**2 + 2*n*m_N**3*math.sqrt(m_N**2 + p_0**2)
    D = b**2 - 4*a*c
    try:
        E_p = (-b/2-math.sqrt(D/4))/a
        #print(p_0, E_p, math.sqrt(p_0**2 + m_N**2) - E_p)
    except:
        E_p = None
    return E_p
def extract_features(imgs, feature_fns, verbose=False):
  """
  Given pixel data for images and several feature functions that can operate on
  single images, apply all feature functions to all images, concatenating the
  feature vectors for each image and storing the features for all images in
  a single matrix.

  Inputs:
  - imgs: N x H X W X C array of pixel data for N images.
  - feature_fns: List of k feature functions. The ith feature function should
    take as input an H x W x D array and return a (one-dimensional) array of
    length F_i.
  - verbose: Boolean; if true, print progress.

  Returns:
  An array of shape (N, F_1 + ... + F_k) where each column is the concatenation
  of all features for a single image.
  """
  num_images = imgs.shape[0]
  if num_images == 0:
    return np.array([])

  # Use the first image to determine feature dimensions
  feature_dims = []
  first_image_features = []
  for feature_fn in feature_fns:
    feats = feature_fn(imgs[0].squeeze())
    assert len(feats.shape) == 1, 'Feature functions must be one-dimensional'
    feature_dims.append(feats.size)
    first_image_features.append(feats)

  # Now that we know the dimensions of the features, we can allocate a single
  # big array to store all features as columns.
  total_feature_dim = sum(feature_dims)
  imgs_features = np.zeros((num_images, total_feature_dim))
  imgs_features[0] = np.hstack(first_image_features).T

  # Extract features for the rest of the images.
  for i in range(1, num_images):
    idx = 0
    for feature_fn, feature_dim in zip(feature_fns, feature_dims):
      next_idx = idx + feature_dim
      imgs_features[i, idx:next_idx] = feature_fn(imgs[i].squeeze())
      idx = next_idx
    if verbose and i % 1000 == 0:
      print('Done extracting features for %d / %d images' % (i, num_images))

  return imgs_features

p_0 = [i/10 for i in range (0,2010)]
E_0 = list(map(energy_finder, p_0))
extract_features(p_0, E_0)
plt.plot(p_0, E_0, label = 'Equation (7)')
plt.xscale('log')
plt.axis((0.1, 200, 0.9, 1.2))
plt.xlabel('$p_0$, GeV')
plt.ylabel('$E_1^*$, GeV')
plt.plot([0, 1000], [1.18, 1.18],'g--', label = '1.18 GeV')
plt.legend()
plt.show()


print(energy_finder(100000000000))
