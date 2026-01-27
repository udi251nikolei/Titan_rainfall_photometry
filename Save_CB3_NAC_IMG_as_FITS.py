import functs
import Get_image_info

image_files = Get_image_info.CB3_NAC_image_files #!!! Need to change camera typr here
image_files_length = len(image_files)
print(f"Total number of image files: {image_files_length}\n")

camera_filter = "CB3"
camera = "NAC" #!!! Need to change camera typr here; always capitalize

total_converted_images = 0
total_failed_conversion = 0
failed_converted_files = []

for n in range(image_files_length):
    print(f"IMG {n+1}: {image_files[n]}")
    try:
        functs.pdr_to_fits_images(image_files[n], camera_filter, camera)
    except Exception as e:
        print(f"{e}\n")
        failed_converted_files.append(str(image_files[n]))
        total_failed_conversion += 1
    else:
        functs.plot_fits_images(image_files[n], camera_filter, camera)
        functs.plot_number_of_rejected_pixels(image_files[n], camera_filter, camera)
        total_converted_images += 1
    
print(f"{camera_filter} {camera} Camera images converted successfully!\n")

print(f"Total number of image files: {image_files_length}")
print(f"Number of successfully converted images: {total_converted_images}")
print(f"Number of failed conversion: {total_failed_conversion}")
print(f"Failed converted files: {failed_converted_files}")

