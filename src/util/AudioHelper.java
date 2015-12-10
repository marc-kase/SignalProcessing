package util;

import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.UnsupportedAudioFileException;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

public class AudioHelper {
        public static void main(String args[]) throws IOException {
            File WAV_FILE = new File("/home/cybersecurity/Desktop/scream2.wav");
            ByteArrayOutputStream out = new ByteArrayOutputStream();
            AudioInputStream in = null;
            try {
                in = AudioSystem.getAudioInputStream(WAV_FILE);
            } catch (UnsupportedAudioFileException e) {

                e.printStackTrace();
            } catch (IOException e) {

                e.printStackTrace();
            }


            int read, i;
            byte[] buff = new byte[1024];
            while ((read = in.read(buff)) > 0) {
                out.write(buff, 0, read);
            }
            out.flush();
            byte[] audioBytes = out.toByteArray();
            System.out.println();

        }
}
