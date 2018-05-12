library(here)
library(hexSticker)

img_path <- here("hex", "vdl-logo-transparent.png")
sticker_path <- here("hex", "tmle3_sticker.png")

sticker(img_path, package = "tmle3",
        p_color = "#242A30", p_size = 5, p_x = 1, p_y = 0.3,
        s_x = 1, s_y = 1.1, s_width = 0.6,
        h_fill = "#ffffff", h_color = "#000000",
        filename = sticker_path)

