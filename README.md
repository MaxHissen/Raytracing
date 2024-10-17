# Raytracing

This is a project I did during my third year at UVic. It was an assignment for my computer graphics class. It was also my first experience programming using c++, and was very helpful in understanding syntax and getting used to commenting on projects and collaborating with peers on getting work done.

The file itself was based on points and lines of the object. Much of the basics of the project was pre-programmed, but there was still a lot to do such as

To get started, the simple shading could be thought of as the shadow of the object. If the rays pass in between the 3 points making a triangle as part of the rasteration, it returns red. Nothing else too it; no shading yet.

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/simple.png?raw=true)

This is the orthographic prespective of the object. Changing to this view style involved starting the camera rays with parallel directions and different positions

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/simple_ortho?raw=true)

Next step was flat shading. This was practice in calculating the plane the triangle exists in, and understanding how rays reflect off of this plane. The brighness of the pixel the ray represents is given by the replection angle, specifically the shallower the reflection, the brighter the reflection.

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/flat_shading.png?raw=true)

This is identical as the previous, but by orthographic perspective.

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/flat_shading_ortho?raw=true)

This animation was mostly pre-programmed for us as part of the project. 

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/flat_shading.gif?raw=true)

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/pv_shading.png?raw=true)
![alt text](https://github.com/MaxHissen/Raytracing/blob/main/pv_shading.gif?raw=true)

![alt text](https://github.com/MaxHissen/Raytracing/blob/main/wireframe.png?raw=true)
![alt text](https://github.com/MaxHissen/Raytracing/blob/main/wireframe.gif?raw=true)
