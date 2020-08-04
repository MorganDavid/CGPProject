#pragma once
template<class M>
class myMatrix {// remeber to set height width and size. 
public:
    M** data = nullptr;
    size_t height = 0;
    size_t width = 0;
};//TODO destructor