palette = {
    'orange': '#e7921e',
    'bluepurple': '#6d78b8',
    'teal': '#99c945',
    'green': '#99c945',
    'pink': '#d170b7',
    'darkgreen': '#39867a',
    'yellow': '#daa51b',
    'blue': '#2f8ac4',
    'purple': '#835fa8',
    'red': '#ed645a',
    'darkpink': '#d14d99',
    'lightgray': '#b7bbad',
    'darkgray': '#84887a',
    'white': '#f6f6f4'
}

def hex_to_rgb(hex):
    hex = hex.lstrip('#')
    rgb = tuple(int(hex[i:i+2], 16) for i in [0, 2, 4])
    return rgb

def pymol_color(hex):
    out = '0x'
    return out + str(hex).lstrip('#')