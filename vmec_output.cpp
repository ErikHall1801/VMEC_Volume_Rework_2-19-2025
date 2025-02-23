

#include "vmec_output.h"


void image_generate(pix_realNumber* pixels_in, pix_RGB* pixels_out)
{
    realNumber color_temp;
    for(int y=0;y<settings.image_height;++y)
    {
        for(int x=0;x<settings.image_width;++x)
        {
            color_temp=std::min(std::max(float(pixels_in[x+y*settings.image_width].R),0.0f),1.0f);
            pixels_out[x+y*settings.image_width].R=static_cast<unsigned char>(255.0f*color_temp);
            color_temp=std::min(std::max(float(pixels_in[x+y*settings.image_width].G),0.0f),1.0f);
            pixels_out[x+y*settings.image_width].G=static_cast<unsigned char>(255.0f*color_temp);
            color_temp=std::min(std::max(float(pixels_in[x+y*settings.image_width].B),0.0f),1.0f);
            pixels_out[x+y*settings.image_width].B=static_cast<unsigned char>(255.0f*color_temp);
        }
    }
}


void write_binary_RBG(pix_RGB* pixels_all, int frame_number) {
    std::string filename = "Render/bin/Render000.bin";
    std::sprintf((char *)filename.data(), "Render/bin/Render%03d.bin", frame_number);
    int img_width=settings.image_width,img_height=settings.image_height;
    std::ofstream fs(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    fs.write(reinterpret_cast<const char*>(&img_width), sizeof(int));
    fs.write(reinterpret_cast<const char*>(&img_height), sizeof(int));
    fs.write(reinterpret_cast<const char*>(pixels_all), settings.image_height*settings.image_width*3);
    fs.close();
}



void write_binary_realNumber(pix_realNumber* pixels_all, int frame_number) {
    std::string filename = "Render/bin/Render000.bin";
    std::sprintf((char *)filename.data(), "Render/bin/Render%03d.bin", frame_number);
    int img_width=settings.image_width,img_height=settings.image_height;
    std::ofstream fs(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    fs.write(reinterpret_cast<const char*>(&img_width), sizeof(int));
    fs.write(reinterpret_cast<const char*>(&img_height), sizeof(int));
    fs.write(reinterpret_cast<const char*>(pixels_all), settings.image_height*settings.image_width*3*sizeof(realNumber));
    fs.close();
}
