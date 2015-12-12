package util;

import javax.sound.sampled.*;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

public class AudioHelper {
    public static AudioFileFormat.Type targetType = AudioFileFormat.Type.SND;

    public static void main(String args[]) throws IOException, UnsupportedAudioFileException {

        String source = "data/Piano.pp.B0.aiff";
        String target = "data/out.txt";

        File file = new File(source);
        AudioInputStream in= AudioSystem.getAudioInputStream(file);
        AudioInputStream ais = AudioSystem.getAudioInputStream(in);

        int read;
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        byte[] buff = new byte[1024];
        while ((read = in.read(buff)) > 0) {
            out.write(buff, 0, read);
        }
        out.flush();

/*        File WAV_FILE = new File("data/Piano.pp.B0.aiff");
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        AudioInputStream in = null;
        try {
            in = AudioSystem.getAudioInputStream(WAV_FILE);
        } catch (UnsupportedAudioFileException | IOException e) {
            e.printStackTrace();
        }


        int read;
        byte[] buff = new byte[1024];
        while ((read = in.read(buff)) > 0) {
            out.write(buff, 0, read);
        }
        out.flush();

        FileWriter fw = new FileWriter(new );
        byte[] audioBytes = out.toByteArray();
        for (int i = 0; i < audioBytes.length; i++) {
            System.out.print((double) audioBytes[i] + ",");
        }*/
    }
}
